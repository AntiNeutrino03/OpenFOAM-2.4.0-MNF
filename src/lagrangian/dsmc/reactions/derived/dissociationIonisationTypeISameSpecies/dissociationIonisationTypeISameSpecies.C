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

Description

\*---------------------------------------------------------------------------*/

#include "dissociationIonisationTypeISameSpecies.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(dissociationIonisationTypeISameSpecies, 0);

addToRunTimeSelectionTable
(
    dsmcReaction,
    dissociationIonisationTypeISameSpecies,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dissociationIonisationTypeISameSpecies::dissociationIonisationTypeISameSpecies
(
    Time& t,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcReaction(t, cloud, dict),
    productIdsIonisation_(),
    productIdsDissociation_(),
    nTotIonisationReactions_(0),
    nIonisationReactionsPerTimeStep_(0),
    nTotDissociationReactions_(0),
    nDissociationReactionsPerTimeStep_(0),
    heatOfReactionIonisation_(),
    heatOfReactionDissociation_()
    {}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dissociationIonisationTypeISameSpecies::
~dissociationIonisationTypeISameSpecies()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dissociationIonisationTypeISameSpecies::initialConfiguration()
{
    setCommonReactionProperties();
    setProperties();
}


void dissociationIonisationTypeISameSpecies::setProperties()
{
    
    if(reactantIds_[0] != reactantIds_[1])
    {
        FatalErrorIn("dissociationIonisationTypeISameSpecies::setProperties()")
            << "Reactants must be same species!"
            << nl
            << exit(FatalError);
    }
    
    if (rDof1_ < VSMALL)
    {
        FatalErrorIn("dissociationIonisationTypeISameSpecies::setProperties()")
            << "Reactants must be molecules: "
            << reactants_[0]
            << nl
            << exit(FatalError);
    }

    if (vDof1_ > 1)
    {
         FatalErrorIn("dissociationIonisationTypeISameSpecies::setProperties()")
            << "Reactions are currently only implemented for "
            << "monatomic and diatomic species."
            << " This is a polyatomic:" << reactants_[0]
            << nl
            << exit(FatalError);
    }

    // reading in ionisation products

    const List<word> productMoleculesIonisation
                    (propsDict_.lookup("productsOfIonisedMolecule"));

    if (productMoleculesIonisation.size() != 2)
    {
        FatalErrorIn("dissociationIonisationTypeISameSpecies::setProperties()")
            << "There should be two products, instead of "
            << productMoleculesIonisation.size() << nl
            << exit(FatalError);
    }

    productIdsIonisation_.setSize(productMoleculesIonisation.size(), -1);

    forAll(productIdsIonisation_, r)
    {
        productIdsIonisation_[r] = findIndex(cloud_.typeIdList(),
                                             productMoleculesIonisation[r]);

        // check that products belong to the typeIdList
        if (productIdsIonisation_[r] == -1)
        {
            FatalErrorIn("dissociationIonisationTypeISameSpecies::setProperties()")
                << "Cannot find type id: " << productMoleculesIonisation[r]
                << nl
                << exit(FatalError);
        }
    }

    scalar rDof2 =
      cloud_.constProps(productIdsIonisation_[0]).rotationalDegreesOfFreedom();

    if (rDof2 < VSMALL)
    {
        FatalErrorIn("dissociationIonisationTypeISameSpecies::setProperties()")
            << "First ionisation product must be a molecular ion: "
            << productMoleculesIonisation[0]
            << nl
            << exit(FatalError);
    }

    const label charge = cloud_.constProps(productIdsIonisation_[1]).charge();

    if (charge != -1)
    {
        FatalErrorIn("dissociationIonisationTypeISameSpecies::setProperties()")
            << "Second ionisation product must be an electron: "
            << productMoleculesIonisation[1]
            << nl
            << exit(FatalError);
    }

    // reading in dissociation products

    const List<word> productMoleculesDissociation
        (propsDict_.lookup("productsOfDissociatedMolecule"));

    if (productMoleculesDissociation.size() != 2)
    {
        FatalErrorIn("dissociationIonisationTypeISameSpecies::setProperties()")
            << "There should be two products, instead of "
            << productMoleculesDissociation.size() << nl
            << exit(FatalError);
    }

    productIdsDissociation_.setSize(productMoleculesDissociation.size(), -1);

    forAll(productIdsDissociation_, r)
    {
        productIdsDissociation_[r] =
        findIndex(cloud_.typeIdList(), productMoleculesDissociation[r]);

        // check that products belong to the typeIdList
        if (productIdsDissociation_[r] == -1)
        {
            FatalErrorIn("dissociationIonisationTypeISameSpecies::setProperties()")
                << "Cannot find type id: " << productMoleculesDissociation[r]
                << nl
                << exit(FatalError);
        }

        // check that products are 'ATOMS' (not 'MOLECULES')

        label iD = productIdsDissociation_[r];

        scalar rDof3 =
                cloud_.constProps(iD).rotationalDegreesOfFreedom();

        if (rDof3 > 1)
        {
            FatalErrorIn("dissociationIonisationTypeISameSpecies::setProperties()")
                << "Dissociation product must be an atom (not a molecule): "
                << productMoleculesDissociation[r]
                << nl
                << exit(FatalError);
        }
    }
    
    heatOfReactionIonisation_ =
        readScalar(propsDict_.lookup("heatOfReactionIonisation"));
        
    heatOfReactionDissociation_ =
        readScalar(propsDict_.lookup("heatOfReactionDissociation"));
}


bool dissociationIonisationTypeISameSpecies::tryReactMolecules
(label typeIdP, const label typeIdQ)
{
    label reactantPId = findIndex(reactantIds_, typeIdP);
    label reactantQId = findIndex(reactantIds_, typeIdQ);

    if (reactantPId == reactantQId)
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


void dissociationIonisationTypeISameSpecies::reaction
(
    dsmcParcel& p,
    dsmcParcel& q,
    const DynamicList<label>& candidateList,
    const List<DynamicList<label> >& candidateSubList,
    const label candidateP,
    const List<label>& whichSubCell
)
{}


void dissociationIonisationTypeISameSpecies::reaction
(
    dsmcParcel& p,
    dsmcParcel& q
)
{
    label typeIdP = p.typeId();
    label typeIdQ = q.typeId();

    if (typeIdP == typeIdQ && typeIdP == reactantIds_[0])
    {
        relax_ = true;
        
        scalar totalReactionProbability = 0.0;
        scalarList reactionProbabilities(4, 0.0);

        bool dissocReactionP = false;
        bool ionisationReactionP = false;
        bool dissocReactionQ = false;
        bool ionisationReactionQ = false;

        // 4 reactions possible
        
        // 1. Dissociation of P
        // 2. Ionisation of P
        // 3. Dissociation of Q
        // 4. Ionisation of Q

        scalar EcPP = 0.0;

        EcPP = reactionProperties(cloud_, p, q).translationalEnergy() + reactionProperties(cloud_, p, q).EVibP();
        label imaxP = EcPP/(physicoChemical::k.value()*reactionProperties(cloud_, p, q).thetaVP());
        label idP = reactionProperties(cloud_, p, q).thetaDP()/reactionProperties(cloud_, p, q).thetaVP();

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
            totalReactionProbability += 1.0;
            reactionProbabilities[1] = 1.0;
        }
        
        scalar EcPQ = reactionProperties(cloud_, p, q).translationalEnergy() + reactionProperties(cloud_, p, q).EVibQ();

        label imaxQ = EcPQ/(physicoChemical::k.value()*reactionProperties(cloud_, p, q).thetaVQ());
        label idQ = reactionProperties(cloud_, p, q).thetaDQ()/reactionProperties(cloud_, p, q).thetaVQ();

        if (imaxQ - idQ > 0)
        {
            // Dissociation can occur
            totalReactionProbability += 1.0;
            reactionProbabilities[2] = 1.0;
        }

        scalar ionisationEnergyQ =
                    cloud_.constProps(typeIdQ).ionisationTemperature()
                    *physicoChemical::k.value();

        // calculate if an ionisation of species Q is possible
        EcPQ = reactionProperties(cloud_, p, q).translationalEnergy() + reactionProperties(cloud_, p, q).EEleQ();

        if ((EcPQ - ionisationEnergyQ) > VSMALL)
        {
            totalReactionProbability += 1.0;
            reactionProbabilities[3] = 1.0;
        }

        // Decide if a reaction is to occur

        if (totalReactionProbability > cloud_.rndGen().scalar01())
        {
            // A chemical reaction is to occur, choose which one

            scalarList normalisedProbabilities(reactionProbabilities.size(),                                                                   0.0);
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
                            dissocReactionP = true;
                            break;
                        }
                        if (i == 1)
                        {
                            // Dissociation reaction is to occur
                            ionisationReactionP = true;
                            break;
                        }
                        if (i == 2)
                        {
                            // Dissociation is to occur
                            dissocReactionQ = true;
                            break;
                        }
                        if (i == 3)
                        {
                            // Ionisation reaction is to occur
                            ionisationReactionQ = true;
                            break;
                        }
                    }
                }
            }
        }

        // dissociation of P
        if (dissocReactionP)
        {
            nTotDissociationReactions_++;
            nDissociationReactionsPerTimeStep_++;

            if (allowSplitting_)
            {
                relax_ = false;
                
                reactionProperties(cloud_, p, q).dissociateP
                (
                    heatOfReactionDissociation_,
                    productIdsDissociation_,
                    p,
                    q
                );
            }
        }

        // dissociation of Q
        if (dissocReactionQ )
        {
            nTotDissociationReactions_++;
            nDissociationReactionsPerTimeStep_++;

            if (allowSplitting_)
            {
                relax_ = false;
                
                reactionProperties(cloud_, p, q).dissociateQ
                (
                    heatOfReactionDissociation_,
                    productIdsDissociation_,
                    p,
                    q
                );
            }
        }

        // ionisation of P
        if (ionisationReactionP)
        {
            nTotIonisationReactions_++;
            nIonisationReactionsPerTimeStep_++;

            if (allowSplitting_)
            {
                relax_ = false;
                
                reactionProperties(cloud_, p, q).ioniseP
                (
                    heatOfReactionIonisation_,
                    productIdsIonisation_,
                    p,
                    q
                );
            }
        }

       // ionisation of Q
        if (ionisationReactionQ)
        {
            nTotIonisationReactions_++;
            nIonisationReactionsPerTimeStep_++;

            if (allowSplitting_)
            {
                relax_ = false;
                
                reactionProperties(cloud_, p, q).ioniseQ
                (
                    heatOfReactionIonisation_,
                    productIdsIonisation_,
                    p,
                    q
                );
            }
        }
    }
}


void  dissociationIonisationTypeISameSpecies::outputResults
(const label counterIndex)
{
    if (writeRatesToTerminal_ == true)
    {
        // measure density

        const List<DynamicList<dsmcParcel*> >& cellOccupancy
            = cloud_.cellOccupancy();

        volume_ = 0.0;

        label molsReactants = 0;

        forAll(cellOccupancy, c)
        {
            const List<dsmcParcel*>& parcelsInCell = cellOccupancy[c];

            forAll(parcelsInCell, pIC)
            {
                dsmcParcel* p = parcelsInCell[pIC];

                if (findIndex(reactantIds_, p->typeId()) != -1)
                {
                    molsReactants++;
                }
            }

            volume_ += mesh_.cellVolumes()[c];
        }

        scalar volume = volume_;
        label nTotReactions = nTotReactions_;

        // Parallel communication
        if (Pstream::parRun())
        {
            reduce(molsReactants, sumOp<label>());
            reduce(volume, sumOp<scalar>());
            reduce(nTotReactions, sumOp<label>());
        }

        numberDensities_[0] = (molsReactants*cloud().nParticle())/volume;
        numberDensities_[1] = (molsReactants*cloud().nParticle())/volume;

        const scalar deltaT = mesh_.time().deltaT().value();

        word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
        word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];

        word productMolA = cloud_.typeIdList()[productIdsDissociation_[0]];
        word productMolB = cloud_.typeIdList()[productIdsDissociation_[1]];

        word productMolC = cloud_.typeIdList()[productIdsIonisation_[0]];
        word productMolD = cloud_.typeIdList()[productIdsIonisation_[1]];

        if ((numberDensities_[0] > 0.0) && (numberDensities_[1] > 0.0))
        {
            scalar reactionRateIonisation = 0.0;
            scalar reactionRateDissociation = 0.0;

            reactionRateIonisation =
            (
                nTotIonisationReactions_
                * cloud_.nParticle()
)/(counterIndex*deltaT*numberDensities_[0]
            * numberDensities_[1]*volume);

            Info<< "Ionisation type I reaction "
                <<  reactantMolA << " + " << reactantMolB
                <<  " --> "
                << productMolC << " + " << productMolD << " + " << reactantMolB
                << ", reaction rate = " << reactionRateIonisation
                << endl;

            reactionRateDissociation =
            (
                nTotDissociationReactions_
                * cloud_.nParticle()
)/(counterIndex*deltaT*numberDensities_[0]
            * numberDensities_[1]*volume);

            Info<< "Dissociation type I reaction "
                <<  reactantMolA << " + " << reactantMolB
                <<  " --> "
                << productMolA << " + " << productMolB << " + " << reactantMolB
                << ", reaction rate = " << reactionRateDissociation
                << endl;
        }
    }
    else
    {
        label nTotDissociationReactions = nTotDissociationReactions_;
        label nTotIonisationReactions = nTotIonisationReactions_;
        label nDissociationReactionsPerTimeStep =
                            nDissociationReactionsPerTimeStep_;
        label nIonisationReactionsPerTimeStep =
                            nIonisationReactionsPerTimeStep_;

        if (Pstream::parRun())
        {
            reduce(nTotDissociationReactions, sumOp<label>());
            reduce(nTotIonisationReactions, sumOp<label>());
            reduce(nDissociationReactionsPerTimeStep, sumOp<label>());
            reduce(nIonisationReactionsPerTimeStep, sumOp<label>());
        }

       if (nTotDissociationReactions > VSMALL)
       {
            word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
            word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];

            word productMolA = cloud_.typeIdList()[productIdsDissociation_[0]];
            word productMolB = cloud_.typeIdList()[productIdsDissociation_[1]];

            Info<< "Dissociation type I reaction "
                <<  reactantMolA << " + " << reactantMolB
                <<  " --> "
                << productMolA << " + " << productMolB << " + " << reactantMolB
                << " is active, nReactions this time step = "
                << nDissociationReactionsPerTimeStep << endl;
        }

        if (nTotIonisationReactions > VSMALL)
        {
            word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
            word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];

            word productMolA = cloud_.typeIdList()[productIdsIonisation_[0]];
            word productMolB = cloud_.typeIdList()[productIdsIonisation_[1]];

            Info<< "Ionisation type I reaction "
                <<  reactantMolA << " + " << reactantMolB
                <<  " --> "
                << productMolA << " + " << productMolB << " + " << reactantMolB
                << " is active, nReactions this time step = "
                << nIonisationReactionsPerTimeStep << endl;
        }
    }

    nDissociationReactionsPerTimeStep_ = 0;
    nIonisationReactionsPerTimeStep_ = 0;
}


bool dissociationIonisationTypeISameSpecies::relax() const
{
    return relax_;
}

} // End namespace Foam

// ************************************************************************* //
