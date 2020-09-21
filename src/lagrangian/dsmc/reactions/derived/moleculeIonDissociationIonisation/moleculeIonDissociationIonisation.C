/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "moleculeIonDissociationIonisation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(moleculeIonDissociationIonisation, 0);

addToRunTimeSelectionTable
(
    dsmcReaction,
    moleculeIonDissociationIonisation,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

moleculeIonDissociationIonisation::moleculeIonDissociationIonisation
(
    Time& t,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcReaction(t, cloud, dict),
    dissociationProducts_(),
    ionisationProducts_(),
    heatOfReactionDissociationAB_(),
    heatOfReactionIonisationAB_(),
    nTotABDissociationReactions_(0),
    nTotABIonisationReactions_(0),
    nDissociationReactionsPerTimeStep_(0),
    nIonisationReationsPerTimeStep_(0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

moleculeIonDissociationIonisation::~moleculeIonDissociationIonisation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void moleculeIonDissociationIonisation::initialConfiguration()
{
    setCommonReactionProperties();
    setProperties();
}


void moleculeIonDissociationIonisation::setProperties()
{
    if (reactantIds_[0] == reactantIds_[1])
    {
        FatalErrorIn("moleculeIonDissociationIonisation::setProperties()")
            << "Reactant molecules cannot be same species." << nl
            << exit(FatalError);
    }

    forAll(reactantIds_, r)
    {
        // check that the first reactant is a 'MOLECULE' (not an 'ATOM')

        if (r == 0)
        {
            label rDof =
               cloud_.constProps(reactantIds_[r]).rotationalDegreesOfFreedom();

            if (rDof < 1)
            {
                FatalErrorIn("moleculeIonDissociationIonisation::setProperties()")
                    << "First reactant must be a molecule (not an atom): "
                    << reactants_[r]
                    << exit(FatalError);
            }

            // check that reactant one only has a single
            // vibrational degree of freedom

            label vDof =
              cloud_.constProps(reactantIds_[r]).vibrationalDegreesOfFreedom();

            if (vDof > 1)
            {
                FatalErrorIn("moleculeIonDissociationIonisation::setProperties()")
                    << "Reactions are currently only implemented "
                    << "for monatomic and diatomic species"
                    << " This is a polyatomic:" << reactants_[r]
                    << nl
                    << exit(FatalError);
            }
        }
    }

    // reading in dissociation products

     List<word> dissociationProducts
     (propsDict_.lookup("productsOfDissociatedMolecule"));

    if (dissociationProducts.size() != 2)
    {
        FatalErrorIn("moleculeIonDissociationIonisation::setProperties()")
            << "number of reactant molecules to be dissociated = "
            << reactantIds_.size()
            << ", should be 2."
            << exit(FatalError);
    }


    dissociationProducts_.setSize(dissociationProducts.size());

    const List<word>& productsForDiss = dissociationProducts;

    if (productsForDiss.size() != 2)
    {
        FatalErrorIn("moleculeIonDissociationIonisation::setProperties()")
            << "There should be two products (for the dissociating molecule "
            << reactants_[0] << "), instead of "
            << productsForDiss.size() << nl
            << exit(FatalError);
    }

    forAll(dissociationProducts_, p)
    {
        dissociationProducts_[p] =
        findIndex(cloud_.typeIdList(), productsForDiss[p]);

        if (dissociationProducts_[p] == -1)
        {
            FatalErrorIn("moleculeIonDissociationIonisation::setProperties()")
                << "Cannot find type id: " << productsForDiss[p] << nl
                << exit(FatalError);
        }
    }

    // reading in ionisation products

    List<word> ionisationProducts
    (propsDict_.lookup("productsOfIonisedMolecule"));

    if (ionisationProducts.size() != 2)
    {
        FatalErrorIn("moleculeIonDissociationIonisation::setProperties()")
            << "number of reactant molecules to be ionised = "
            << reactantIds_.size()
            << ", should be 2."
            << exit(FatalError);
    }


    ionisationProducts_.setSize(ionisationProducts.size());

    const List<word>& productsForIon = ionisationProducts;

    if (productsForIon.size() != 2)
    {
        FatalErrorIn("moleculeIonDissociationIonisation::setProperties()")
            << "There should be two products (for the ionising molecule "
            << reactants_[0] << "), instead of "
            << productsForIon.size() << nl
            << exit(FatalError);
    }

    ionisationProducts_.setSize(productsForIon.size(), -1);

    forAll(ionisationProducts_, p)
    {
        ionisationProducts_[p] =
        findIndex(cloud_.typeIdList(), productsForIon[p]);

        if (ionisationProducts_[p] == -1)
        {
            FatalErrorIn("moleculeIonDissociationIonisation::setProperties()")
                << "Cannot find type id: " << productsForIon[p] << nl
                << exit(FatalError);
        }
    }

    // check ionisation product 2 is an electron
    forAll(ionisationProducts_, p)
    {
        if (p == 1)
        {
            const label charge =
            cloud_.constProps(ionisationProducts_[p]).charge();

            if (charge != -1)
            {
                FatalErrorIn("moleculeIonDissociationIonisation::setProperties()")
                    << "Second ionisation product must be an electron: "
                    << productsForIon[p] << nl
                    << exit(FatalError);
            }
        }
    }
    
    heatOfReactionDissociationAB_ =
        readScalar(propsDict_.lookup("heatOfReactionDissociation"));
    heatOfReactionIonisationAB_ =
        readScalar(propsDict_.lookup("heatOfReactionIonisation"));
}


bool moleculeIonDissociationIonisation::
tryReactMolecules(label typeIdP, const label typeIdQ)
{
    label reactantPId = findIndex(reactantIds_, typeIdP);
    label reactantQId = findIndex(reactantIds_, typeIdQ);

    if (reactantPId != reactantQId)
    {
        if
        (
            (reactantPId != -1) &&
            (reactantQId != -1)
        )
        {
            return true;
        }
    }

    if
    (
        (reactantPId != -1) &&
        (reactantQId != -1) &&
        (reactantPId != reactantQId)
    )
    {
        return true;
    }
    else
    {
        return false;
    }
}


void moleculeIonDissociationIonisation::reaction
(
    dsmcParcel& p,
    dsmcParcel& q,
    const DynamicList<label>& candidateList,
    const List<DynamicList<label> >& candidateSubList,
    const label candidateP,
    const List<label>& whichSubCell
)
{}


void moleculeIonDissociationIonisation::reaction
(
    dsmcParcel& p,
    dsmcParcel& q
)

{
    label typeIdP = p.typeId();
    label typeIdQ = q.typeId();

    if (typeIdP == reactantIds_[0] && typeIdQ == reactantIds_[1])
    {
        scalar totalReactionProbability = 0.0;
        scalarList reactionProbabilities(2, 0.0);

        relax_ = true;

        bool dissocReactionP = false;
        bool ionisationReactionP = false;

        // 2 reactions possible
        
        // 1. Dissociation of P
        // 2. Ionisation of P

        scalar EcPP = 0.0;
        label idP = cloud_.constProps(typeIdP).charDissQuantumLevel()[0];
        label imaxP = 0;

        // calculate if a dissociation of species P is possible
        EcPP = reactionProperties(cloud_, p, q).translationalEnergy() + reactionProperties(cloud_, p, q).EVibP();

        imaxP = EcPP/(physicoChemical::k.value()*reactionProperties(cloud_, p, q).thetaVP());

        if (imaxP - idP > 0)
        {
            // Dissociation can occur
            totalReactionProbability += 1.0;
            reactionProbabilities[0] = 1.0;
        }

        scalar ionisationEnergy =
                    cloud_.constProps(typeIdP).ionisationTemperature()
                    *physicoChemical::k.value();

        // calculate if an ionisation of species P is possible
        EcPP = reactionProperties(cloud_, p, q).translationalEnergy() + reactionProperties(cloud_, p, q).EEleP();

        if ((EcPP - ionisationEnergy) > VSMALL)
        {
            // Ionisation can occur
//             totalReactionProbability += 1.0;
//             reactionProbabilities[1] = 1.0;
        }

        // Decide if a reaction is to occur

        if (totalReactionProbability > cloud_.rndGen().scalar01())
        {
            // A chemical reaction is to occur, choose which one

            scalarList normalisedProbabilities(reactionProbabilities.size(),
                                                            0.0);
            scalar cumulativeProbability = 0.0;

            normalisedProbabilities = reactionProbabilities
                    /totalReactionProbability;

            forAll(normalisedProbabilities, i)
            {
                // If current reaction can't occur, don't check for it
                if (normalisedProbabilities[i] > VSMALL)
                {
                    cumulativeProbability += normalisedProbabilities[i];

                    if (cumulativeProbability > cloud_.rndGen().scalar01())
                    {
                        // Current reaction is to occur

                        if (i == 0)
                        {
                            // Dissociation is to occur
                            dissocReactionP = true;
                            break;
                        }
                        if (i == 1)
                        {
                            // Ionisation reaction is to occur
                            ionisationReactionP = true;
                            break;
                        }
                    }
                }
            }
        }

        if (dissocReactionP)
        {
            nDissociationReactionsPerTimeStep_++;
            nTotABDissociationReactions_++;

            if (allowSplitting_)
            {
                relax_ = false;
                
                reactionProperties(cloud_, p, q).dissociateP
                (
                    heatOfReactionDissociationAB_,
                    dissociationProducts_,
                    p,
                    q
                );
            }
        }

        if (ionisationReactionP)
        {
            nIonisationReationsPerTimeStep_++;
            nTotABIonisationReactions_++;

            if (allowSplitting_)
            {
                relax_ = false;
                
                reactionProperties(cloud_, p, q).ioniseP
                (
                    heatOfReactionIonisationAB_,
                    ionisationProducts_,
                    p,
                    q
                );
            }
        }
    }

    if (typeIdP == reactantIds_[1] && typeIdQ == reactantIds_[0])
    {
        relax_ = true;

        scalar totalReactionProbability = 0.0;
        scalarList reactionProbabilities(2, 0.0);

        bool dissocReactionQ = false;
        bool ionisationReactionQ = false;

        // 2 reactions possible
        
        // 1. Dissociation of Q
        // 2. Ionisation of Q

        scalar EcQP = 0.0;
        label idQ = cloud_.constProps(typeIdQ).charDissQuantumLevel()[0];
        label imaxQ = 0;

        // calculate if a dissociation of species Q is possible
        EcQP = reactionProperties(cloud_, p, q).translationalEnergy() + reactionProperties(cloud_, p, q).EVibQ();

        imaxQ = EcQP/(physicoChemical::k.value()*reactionProperties(cloud_, p, q).thetaVQ());

        if (imaxQ - idQ > 0)
        {
            // Dissociation can occur
            totalReactionProbability += 1.0;
            reactionProbabilities[0] = 1.0;
        }

        scalar ionisationEnergy =
                    cloud_.constProps(typeIdQ).ionisationTemperature()
                    *physicoChemical::k.value();

        // calculate if an ionisation of species Q is possible
        EcQP = reactionProperties(cloud_, p, q).translationalEnergy() + reactionProperties(cloud_, p, q).EEleQ();

        if ((EcQP - ionisationEnergy) > VSMALL)
        {
//             totalReactionProbability += 1.0;
//             reactionProbabilities[1] = 1.0;
        }
        // Decide if a reaction is to occur

        if (totalReactionProbability > cloud_.rndGen().scalar01())
        {
            // A chemical reaction is to occur, choose which one

            scalarList normalisedProbabilities(reactionProbabilities.size(),
                                               0.0);
            scalar cumulativeProbability = 0.0;

            normalisedProbabilities = reactionProbabilities
                                        /totalReactionProbability;

            forAll(normalisedProbabilities, i)
            {
                // If current reaction can't occur, don't check for it
                if (normalisedProbabilities[i] > VSMALL)
                {
                    cumulativeProbability += normalisedProbabilities[i];

                    if (cumulativeProbability > cloud_.rndGen().scalar01())
                    {
                        // Current reaction is to occur

                        if (i == 0)
                        {
                            // Ionisation is to occur
                            dissocReactionQ = true;
                            break;
                        }
                        if (i == 1)
                        {
                            // Dissociation reaction is to occur
                            ionisationReactionQ = true;
                            break;
                        }
                    }
                }
            }
        }

        if (dissocReactionQ)
        {
            nDissociationReactionsPerTimeStep_++;
            nTotABDissociationReactions_++;

            if (allowSplitting_)
            {
                relax_ = false;
                
                reactionProperties(cloud_, p, q).dissociateQ
                (
                    heatOfReactionDissociationAB_,
                    dissociationProducts_,
                    p,
                    q
                );
            }
        }

        if (ionisationReactionQ)
        {
            nIonisationReationsPerTimeStep_++;
            nTotABIonisationReactions_++;

            if (allowSplitting_)
            {
                relax_ = false;
                
                reactionProperties(cloud_, p, q).ioniseQ
                (
                    heatOfReactionIonisationAB_,
                    ionisationProducts_,
                    p,
                    q
                );
            }
        }
    }
}


void  moleculeIonDissociationIonisation::
outputResults(const label counterIndex)
{
    if (writeRatesToTerminal_ == true)
    {
        // measure density
        const List<DynamicList<dsmcParcel*> >& cellOccupancy
            = cloud_.cellOccupancy();

        volume_ = 0.0;

        List<label> mols (2, 0);

        forAll(cellOccupancy, c)
        {
            const List<dsmcParcel*>& parcelsInCell = cellOccupancy[c];

            forAll(parcelsInCell, pIC)
            {
                dsmcParcel* p = parcelsInCell[pIC];

                label id = findIndex(reactantIds_, p->typeId());

                if (id != -1)
                {
                    mols[id]++;
                }
            }

            volume_ += mesh_.cellVolumes()[c];
        }

        scalar volume = volume_;
        label nTotABDissociationReactions = nTotABDissociationReactions_;
        label nTotABIonisationReactions = nTotABIonisationReactions_;

        // Parallel communication
        if (Pstream::parRun())
        {
            reduce(volume, sumOp<scalar>());
            reduce(mols[0], sumOp<label>());
            reduce(mols[1], sumOp<label>());
            reduce(nTotABDissociationReactions, sumOp<label>());
            reduce(nTotABIonisationReactions, sumOp<label>());
        }

        numberDensities_[0] = (mols[0]*cloud().nParticle())/volume;
        numberDensities_[1] = (mols[1]*cloud().nParticle())/volume;

        const scalar deltaT = mesh_.time().deltaT().value();

        word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
        word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];

        word dissociationProductMolA =
            cloud_.typeIdList()[dissociationProducts_[0]];
        word dissociationProductMolB =
            cloud_.typeIdList()[dissociationProducts_[1]];

        word ionisationProductMolA =
            cloud_.typeIdList()[ionisationProducts_[0]];
        word ionisationProductMolB =
            cloud_.typeIdList()[ionisationProducts_[1]];

        if ((numberDensities_[0] > 0.0) && (numberDensities_[1] > 0.0))
        {
            scalar reactionRate1 = 0.0;
            scalar reactionRate2 = 0.0;

            reactionRate1 =
            (
                nTotABDissociationReactions
                * cloud_.nParticle()
)/(counterIndex*deltaT*numberDensities_[0]
            * numberDensities_[1]*volume);

            reactionRate2 =
            (
                nTotABIonisationReactions
                * cloud_.nParticle()
)/(counterIndex*deltaT*numberDensities_[0]
            * numberDensities_[1]*volume);

            Info << "Dissociation type I reaction "
            <<  reactantMolA
            << " + " << reactantMolB << " --> "
            << dissociationProductMolA << " + "
            << dissociationProductMolB << " + "
            << reactantMolB << ", reaction rate = "
            << reactionRate1  << nl
            << "Ionisation type I reaction "
            <<  reactantMolA << " + "
            << reactantMolB << " --> "
            << ionisationProductMolA << " + "
            << ionisationProductMolB << " + "
            << reactantMolB << ", reaction rate = "
            << reactionRate2
            << endl;
        }
    }
    else
    {
        label nTotABDissociationReactions = nTotABDissociationReactions_;
        label nTotABIonisationReactions = nTotABIonisationReactions_;

        label nDissociationReactionsPerTimeStep =
                        nDissociationReactionsPerTimeStep_;
        label nIonisationReationsPerTimeStep =
                            nIonisationReationsPerTimeStep_;

        // Parallel communication
        if (Pstream::parRun())
        {
            reduce(nTotABDissociationReactions, sumOp<label>());
            reduce(nTotABIonisationReactions, sumOp<label>());
            reduce(nDissociationReactionsPerTimeStep, sumOp<label>());
            reduce(nIonisationReationsPerTimeStep, sumOp<label>());
        }

        if (nTotABDissociationReactions > VSMALL)
        {
            word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
            word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];

            word dissociationProductMolA =
                    cloud_.typeIdList()[dissociationProducts_[0]];
            word dissociationProductMolB =
                cloud_.typeIdList()[dissociationProducts_[1]];

            Info << "Dissociation type I reaction " <<  reactantMolA
            << " + " << reactantMolB << " --> "
            << dissociationProductMolA << " + "
            << dissociationProductMolB << " + "
            << reactantMolB
            << " is active, nReactions this time step = "
            << nDissociationReactionsPerTimeStep << endl;
        }

        if (nTotABIonisationReactions > VSMALL)
        {
            word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
            word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];

            word ionisationProductMolA =
                cloud_.typeIdList()[ionisationProducts_[0]];
            word ionisationProductMolB =
                cloud_.typeIdList()[ionisationProducts_[1]];

            Info << "Ionisation type I reaction " <<  reactantMolA
            << " + " << reactantMolB << " --> "
            << ionisationProductMolA << " + "
            << ionisationProductMolB << " + "
            << reactantMolB
            << " is active, nReactions this time step = "
            << nIonisationReationsPerTimeStep << endl;
        }

    }

    nDissociationReactionsPerTimeStep_ = 0.0;
    nIonisationReationsPerTimeStep_ = 0.0;
}


bool moleculeIonDissociationIonisation::relax() const
{
    return relax_;
}

} // End namespace Foam

// ************************************************************************* //
