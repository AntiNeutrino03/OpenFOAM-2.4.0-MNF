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

#include "moleculeElectronDissociationIonisation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(moleculeElectronDissociationIonisation, 0);

addToRunTimeSelectionTable
(
    dsmcReaction,
    moleculeElectronDissociationIonisation,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
moleculeElectronDissociationIonisation::moleculeElectronDissociationIonisation
(
    Time& t,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcReaction(t, cloud, dict),
    productIdsDiss_(),
    productIdsIon_(),
    heatOfReactionDiss_(),
    heatOfReactionIon_(),
    nDissociationReactions_(0),
    nDissociationReactionsPerTimeStep_(0),
    nIonisationReactions_(0),
    nIonisationReactionsPerTimeStep_(0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

moleculeElectronDissociationIonisation::
~moleculeElectronDissociationIonisation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void moleculeElectronDissociationIonisation::initialConfiguration()
{
    setCommonReactionProperties();
    setProperties();
}


void moleculeElectronDissociationIonisation::setProperties()
{
    if (reactantIds_[0] == reactantIds_[1])
    {
        FatalErrorIn("moleculeElectronDissociationIonisation::setProperties()")
            << "Reactant molecules cannot be same species." << nl
            << exit(FatalError);
    }

    // check that reactant one is a 'MOLECULE''

    label rDof1 =
        cloud_.constProps(reactantIds_[0]).rotationalDegreesOfFreedom();

    if (rDof1 < 1)
    {
        FatalErrorIn("moleculeElectronDissociationIonisation::setProperties()")
            << "First reactant must be a molecule "
            << "(not an atom or an electron): " << reactants_[0]
            << nl
            << exit(FatalError);
    }

    label vDof =
            cloud_.constProps(reactantIds_[0]).vibrationalDegreesOfFreedom();

    if (vDof > 1)
    {
        FatalErrorIn("moleculeElectronDissociationIonisation::setProperties()")
            << "Reactions are currently only implemented for "
            << "monatomic and diatomic species"
            << " This is a polyatomic:" << reactants_[0]
            << nl
            << exit(FatalError);
    }

    // check that reactant two is an 'ELECTRON'

    label charge = cloud_.constProps(reactantIds_[1]).charge();

    if (charge != -1)
    {
        FatalErrorIn("moleculeElectronDissociationIonisation::setProperties()")
            << "Second reactant must be an electron, not "
            << reactants_[1]
            << nl
            << exit(FatalError);
    }

    // reading in dissociation products

    List<word> productMoleculesDissociation
    (propsDict_.lookup("productsOfDissociatedMolecule"));

    if (productMoleculesDissociation.size() != 2)
    {
        FatalErrorIn("moleculeElectronDissociationIonisation::setProperties()")
            << "Number of dissociation products is "
            << productMoleculesDissociation.size()
            << ", should be two."
            << exit(FatalError);
    }


    productIdsDiss_.setSize(productMoleculesDissociation.size());

    forAll(productMoleculesDissociation, r)
    {
        if (productIdsDiss_.size() != 2)
        {
            FatalErrorIn("moleculeElectronDissociationIonisation::setProperties()")
                << "There should be two products "
                << "(for the dissociating molecule "
                << reactants_[r] << "), instead of "
                << productIdsDiss_.size() << nl
                << exit(FatalError);
        }

        forAll(productIdsDiss_, r)
        {
            productIdsDiss_[r] =
            findIndex(cloud_.typeIdList(), productMoleculesDissociation[r]);

            // check that reactants belong to the typeIdList
            if (productIdsDiss_[r] == -1)
            {
                FatalErrorIn("moleculeElectronDissociationIonisation::setProperties()")
                    << "Cannot find type id: "
                    << productMoleculesDissociation[r] << nl
                    << exit(FatalError);
            }
        }

        // check that product one is an 'ATOM' (not a 'MOLECULE')

        scalar rDof3 =
            cloud_.constProps(productIdsDiss_[0]).rotationalDegreesOfFreedom();

        if (rDof3 != 0)
        {
            FatalErrorIn("moleculeElectronDissociationIonisation::setProperties()")
                << "First product must be an atom (not a molecule): "
                << productMoleculesDissociation[0]
                << nl
                << exit(FatalError);
        }

        // check that product two is an 'ATOM' (not a 'MOLECULE')

        scalar rDof4 =
            cloud_.constProps(productIdsDiss_[1]).rotationalDegreesOfFreedom();

        if (rDof4 != 0)
        {
            FatalErrorIn("moleculeElectronDissociationIonisation::setProperties()")
                << "Second product must be an atom (not a molecule): "
                << productMoleculesDissociation[1]
                << nl
                << exit(FatalError);
        }
    }

    // reading in ionisation products

    List<word> productMoleculesIonisation
    (propsDict_.lookup("productsOfIonisedMolecule"));

    if (productMoleculesIonisation.size() != 2)
    {
        FatalErrorIn("moleculeElectronDissociationIonisation::setProperties()")
            << "Number of dissociation products is "
            << productMoleculesIonisation.size()
            << ", should be two."
            << exit(FatalError);
    }


    productIdsIon_.setSize(productMoleculesIonisation.size());

    forAll(productMoleculesIonisation, r)
    {
        if (productIdsIon_.size() != 2)
        {
            FatalErrorIn("moleculeElectronDissociationIonisation::setProperties()")
                << "There should be two products "
                << "(for the dissociating molecule "
                << reactants_[r] << "), instead of "
                << productIdsIon_.size() << nl
                << exit(FatalError);
        }

        forAll(productIdsIon_, r)
        {
            productIdsIon_[r] =
            findIndex(cloud_.typeIdList(), productMoleculesIonisation[r]);

            // check that reactants belong to the typeIdList
            if (productIdsIon_[r] == -1)
            {
                FatalErrorIn("moleculeElectronDissociationIonisation::setProperties()")
                    << "Cannot find type id: "
                    << productMoleculesIonisation[r] << nl
                    << exit(FatalError);
            }
        }

        // check that product one is a 'MOLECULE', not an 'ATOM'

        scalar rDof5 =
            cloud_.constProps(productIdsIon_[0]).rotationalDegreesOfFreedom();

        if (rDof5 < 1)
        {
            FatalErrorIn("moleculeElectronDissociationIonisation::setProperties()")
                << "First product must be a molecule ("
                << "not an atom): "
                << productMoleculesIonisation[0]
                << nl
                << exit(FatalError);
        }

        // check that product two is an 'ELECTRON'

        const label charge = cloud_.constProps(productIdsIon_[1]).charge();

        if (charge != -1)
        {
            FatalErrorIn("moleculeElectronDissociationIonisation::setProperties()")
                << "Second product must be an electron: "
                << productMoleculesIonisation[1]
                << nl
                << exit(FatalError);
        }
    }
    
    heatOfReactionDiss_ =
        readScalar(propsDict_.lookup("heatOfReactionDissociation"));
    heatOfReactionIon_ =
        readScalar(propsDict_.lookup("heatOfReactionIonisation"));
}


bool moleculeElectronDissociationIonisation::
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


void moleculeElectronDissociationIonisation::reaction
(
    dsmcParcel& p,
    dsmcParcel& q,
    const DynamicList<label>& candidateList,
    const List<DynamicList<label> >& candidateSubList,
    const label candidateP,
    const List<label>& whichSubCell
)
{}


void moleculeElectronDissociationIonisation::reaction
(
    dsmcParcel& p,
    dsmcParcel& q
)
{
    label typeIdP = p.typeId();
    label typeIdQ = q.typeId();

    if (typeIdP == reactantIds_[0] && typeIdQ == reactantIds_[1])
    {
        relax_ = true;

        scalar totalReactionProbability = 0.0;
        scalarList reactionProbabilities(2, 0.0);
        
        bool dissocReaction = false;
        bool ionisationReaction = false;

        scalar EcPPIon = 0.0;
        scalar ionisationEnergy =
                cloud_.constProps(typeIdP).ionisationTemperature()
                *physicoChemical::k.value();

        // calculate if an ionisation of species P is possible
        EcPPIon = reactionProperties(cloud_, p, q).translationalEnergy() + reactionProperties(cloud_, p, q).EEleP();

        if ((EcPPIon - ionisationEnergy) > VSMALL)
        {
            totalReactionProbability += 1.0;
            reactionProbabilities[0] = 1.0;
        }

        scalar EcPPDiss = 0.0;
        label idP = reactionProperties(cloud_, p, q).charDissLevelP();
        label imaxP = 0;

        // calculate if a dissociation of species P is possible
        EcPPDiss = reactionProperties(cloud_, p, q).translationalEnergy() + reactionProperties(cloud_, p, q).EVibP();

        imaxP = EcPPDiss/(physicoChemical::k.value()*reactionProperties(cloud_, p, q).thetaVP());

        if (imaxP - idP > 0)
        {
            totalReactionProbability += 1.0;
            reactionProbabilities[1] = 1.0;
        }

        // Decide if a reaction is to occur

        if (totalReactionProbability > cloud_.rndGen().scalar01())
        {
            // A chemical reaction is to occur, choose which one

            scalarList normalisedProbabilities(reactionProbabilities.size(),
                                               0.0);
            scalar cumulativeProbability = 0.0;

            forAll(normalisedProbabilities, i)
            {
                normalisedProbabilities[i] = reactionProbabilities[i]
                                                /totalReactionProbability;

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
                            ionisationReaction = true;
                            break;
                        }
                        if (i == 1)
                        {
                            // Dissociation reaction is to occur
                            dissocReaction = true;
                            break;
                        }
                    }
                }
            }
        }

        // Perform a dissociation reaction
        if (dissocReaction)
        {
            nDissociationReactions_++;
            nDissociationReactionsPerTimeStep_++;

            if (allowSplitting_)
            {
                relax_ = false;
                
                reactionProperties(cloud_, p, q).dissociateP
                (
                    heatOfReactionDiss_,
                    productIdsDiss_,
                    p,
                    q
                );
            }
        }

        // Perform an ionisation reaction

        if (ionisationReaction)
        {
            nIonisationReactions_++;
            nIonisationReactionsPerTimeStep_++;

            if (allowSplitting_)
            {
                relax_ = false;
                
                reactionProperties(cloud_, p, q).ioniseP
                (
                    heatOfReactionIon_,
                    productIdsIon_,
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

        bool dissocReaction = false;
        bool ionisationReaction = false;

        scalar EcPQIon = 0.0;
        scalar ionisationEnergy =
                cloud_.constProps(typeIdQ).ionisationTemperature()
                *physicoChemical::k.value();

        // calculate if an ionisation of species Q is possible
        EcPQIon = reactionProperties(cloud_, p, q).translationalEnergy() + reactionProperties(cloud_, p, q).EEleQ();

        if ((EcPQIon - ionisationEnergy) > VSMALL)
        {
            totalReactionProbability += 1.0;
            reactionProbabilities[0] = 1.0;
        }

        scalar EcPQ = 0.0;
        label idQ = reactionProperties(cloud_, p, q).charDissLevelQ();
        label imaxQ = 0;

        // calculate if a dissociation of species Q is possible
        EcPQ = reactionProperties(cloud_, p, q).translationalEnergy() + reactionProperties(cloud_, p, q).EVibQ();

        imaxQ = EcPQ/(physicoChemical::k.value()*reactionProperties(cloud_, p, q).thetaVQ());

        if (imaxQ - idQ > 0)
        {
            totalReactionProbability += 1.0;
            reactionProbabilities[1] = 1.0;
        }

        // Decide if a reaction is to occur

        if (totalReactionProbability > cloud_.rndGen().scalar01())
        {
            // A chemical reaction is to occur, choose which one

            scalarList normalisedProbabilities(reactionProbabilities.size(),
                                               0.0);
            scalar cumulativeProbability = 0.0;

            forAll(normalisedProbabilities, i)
            {
                normalisedProbabilities[i] = reactionProbabilities[i]
                                                /totalReactionProbability;

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
                            ionisationReaction = true;
                            break;
                        }
                        if (i == 1)
                        {
                            // Dissociation reaction is to occur
                            dissocReaction = true;
                            break;
                        }
                    }
                }
            }
        }

        // Perform a dissociation reaction
        if (dissocReaction)
        {
            nDissociationReactions_++;
            nDissociationReactionsPerTimeStep_++;

            if (allowSplitting_)
            {
                relax_ = false;
                
                reactionProperties(cloud_, p, q).dissociateQ
                (
                    heatOfReactionDiss_,
                    productIdsDiss_,
                    p,
                    q
                );
            }
        }

        // Perform an ionisation reaction
        if (ionisationReaction)
        {
            nIonisationReactions_++;
            nIonisationReactionsPerTimeStep_++;

            if (allowSplitting_)
            {
                relax_ = false;
                
                reactionProperties(cloud_, p, q).ioniseQ
                (
                    heatOfReactionIon_,
                    productIdsIon_,
                    p,
                    q
                );
            }
        }
    }
}


void  moleculeElectronDissociationIonisation::
outputResults(const label counterIndex)
{
    if (writeRatesToTerminal_ == true)
    {
        const List<DynamicList<dsmcParcel*> >& cellOccupancy
            = cloud_.cellOccupancy();

        List<label> mols (2, 0);
        volume_ = 0.0;

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
        label nTotReactionsDissociation = nDissociationReactions_;
        label nTotReactionsIonisation = nIonisationReactions_;

        // Parallel communication
        if (Pstream::parRun())
        {
            reduce(volume, sumOp<scalar>());
            reduce(mols[0], sumOp<label>());
            reduce(mols[1], sumOp<label>());
            reduce(nTotReactionsDissociation, sumOp<label>());
            reduce(nTotReactionsIonisation, sumOp<label>());
        }

        numberDensities_[0] = (mols[0]*cloud().nParticle())/volume;
        numberDensities_[1] = (mols[1]*cloud().nParticle())/volume;

        const scalar deltaT = mesh_.time().deltaT().value();

        word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
        word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];

        word productMolA = cloud_.typeIdList()[productIdsDiss_[0]];
        word productMolB = cloud_.typeIdList()[productIdsDiss_[1]];

        word productMolC = cloud_.typeIdList()[productIdsIon_[0]];
        word productMolD = cloud_.typeIdList()[productIdsIon_[1]];

        if ((numberDensities_[0] > 0.0) && (numberDensities_[1] > 0.0))
        {
            scalar reactionRateDissociation = 0.0;

            reactionRateDissociation =
            (
            nTotReactionsDissociation
            * cloud_.nParticle()
)/(counterIndex*deltaT*numberDensities_[0]
            * numberDensities_[1]*volume);

            Info<< "Electron dissociation reaction "
                <<  reactantMolA << " + " << reactantMolB
                <<  " --> "
                << productMolA << " + " << productMolB << " + " << reactantMolB
                << ", reaction rate = " << reactionRateDissociation
                << endl;

            scalar reactionRateIonisation = 0.0;

            reactionRateIonisation =
            (
            nTotReactionsIonisation
            * cloud_.nParticle()
)/(counterIndex*deltaT*numberDensities_[0]
            * numberDensities_[1]*volume);

            Info<< "Electron ionisation reaction "
                <<  reactantMolA << " + " << reactantMolB
                <<  " --> "
                << productMolC << " + " << productMolD << " + " << reactantMolB
                << ", reaction rate = " << reactionRateIonisation
                << endl;
        }
    }
    else
    {
        label nTotReactionsDissociation = nDissociationReactions_;
        label nTotReactionsIonisation = nIonisationReactions_;
        label nDissociationReactionsPerTimeStep =
            nDissociationReactionsPerTimeStep_;
        label nIonisationReactionsPerTimeStep =
            nIonisationReactionsPerTimeStep_;

        if (Pstream::parRun())
        {
            reduce(nTotReactionsDissociation, sumOp<label>());
            reduce(nTotReactionsIonisation, sumOp<label>());

            reduce(nDissociationReactionsPerTimeStep, sumOp<label>());
            reduce(nIonisationReactionsPerTimeStep, sumOp<label>());
        }

        if (nTotReactionsDissociation > VSMALL)
        {
                word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
                word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];

                word productMolA = cloud_.typeIdList()[productIdsDiss_[0]];
                word productMolB = cloud_.typeIdList()[productIdsDiss_[1]];

                Info<< "Electron dissociation reaction "
                <<  reactantMolA << " + " << reactantMolB
                <<  " --> "
                << productMolA << " + " << productMolB << " + "
                << reactantMolB
                << " is active, nReactions this time step = "
                << nDissociationReactionsPerTimeStep << endl;
        }

        if (nTotReactionsIonisation > VSMALL)
        {
                word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
                word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];

                word productMolA = cloud_.typeIdList()[productIdsIon_[0]];
                word productMolB = cloud_.typeIdList()[productIdsIon_[1]];

                Info<< "Electron ionisation reaction "
                <<  reactantMolA << " + " << reactantMolB
                <<  " --> "
                << productMolA << " + " << productMolB << " + "
                << reactantMolB
                << " is active, nReactions this time step = "
                << nIonisationReactionsPerTimeStep << endl;
        }
    }

    nReactionsPerTimeStep_ = 0.0;
    nDissociationReactionsPerTimeStep_ = 0.0;
    nIonisationReactionsPerTimeStep_ = 0.0;
}


bool moleculeElectronDissociationIonisation::relax() const
{
    return relax_;
}

} // End namespace Foam

// ************************************************************************* //
