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

#include "moleculeAtomDissociationIonisation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(moleculeAtomDissociationIonisation, 0);

addToRunTimeSelectionTable
(
    dsmcReaction,
    moleculeAtomDissociationIonisation,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
moleculeAtomDissociationIonisation::moleculeAtomDissociationIonisation
(
    Time& t,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcReaction(t, cloud, dict),
    productIdsDiss_(),
    productIdsIon_(),
    productIdsIon2_(),
    heatOfReactionDiss_(),
    heatOfReactionIon_(),
    heatOfReactionIon2_(),
    nTotReactionsDiss_(0),
    nTotReactionsIon_(0),
    nTotReactionsIon2_(0),
    nDissociationReactionsPerTimeStep_(0),
    nIonisationReactionsPerTimeStep_(0),
    nIonisationReactions2PerTimeStep_(0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

moleculeAtomDissociationIonisation::~moleculeAtomDissociationIonisation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void moleculeAtomDissociationIonisation::initialConfiguration()
{
    setCommonReactionProperties();
    setProperties();
}


void moleculeAtomDissociationIonisation::setProperties()
{
    if (reactantIds_[0] == reactantIds_[1])
    {
        FatalErrorIn("moleculeAtomDissociationIonisation::setProperties()")
            << "Reactant molecules cannot be same species." << nl
            << exit(FatalError);
    }

    // check that first reactant is a 'MOLECULE' (not an 'ATOM')

    label rDof1 =
    cloud_.constProps(reactantIds_[0]).rotationalDegreesOfFreedom();
    
    if (rDof1 < 1)
    {
        FatalErrorIn("moleculeAtomDissociationIonisation::setProperties()")
            << "First reactant must be a molecule (not an atom): "
            << reactants_[0]
            << nl
            << exit(FatalError);
    }
    
    label vDof1 =
    cloud_.constProps(reactantIds_[0]).vibrationalDegreesOfFreedom();

    if (vDof1 > 1)
    {
        FatalErrorIn("moleculeAtomDissociationIonisation::setProperties()")
            << "Reactions are currently only implemented for "
            << "monatomic and diatomic species"
            << " This is a polyatomic:" << reactantIds_[0]
            << nl
            << exit(FatalError);
    }

    // check that second reactant is an 'ATOM' (not a 'MOLECULE')
    
    label rDof2 =
        cloud_.constProps(reactantIds_[1]).rotationalDegreesOfFreedom();
        
    if (rDof2 > VSMALL)
    {
        FatalErrorIn("moleculeAtomDissociationIonisation::setProperties()")
            << "Second reactant must be atom (not a molecule): "
            << reactants_[1]
            << nl
            << exit(FatalError);
    }
    
    // reading in products

    List<word> productMoleculesDiss
    (
        propsDict_.lookup("productsOfDissociatedMolecule")
    );

    if (productMoleculesDiss.size() != 2)
    {
        FatalErrorIn("moleculeAtomDissociationIonisation::setProperties()")
            << "Number of dissociation products is "
            << productMoleculesDiss.size()
            << ", should be two."
            << exit(FatalError);
    }

    List<word> productMoleculesIon
    (
        propsDict_.lookup("productsOfIonisedMolecule")
    );

    if (productMoleculesIon.size() != 2)
    {
        FatalErrorIn("moleculeAtomDissociationIonisation::setProperties()")
            << "Number of ionisation products is "
            << productMoleculesIon.size()
            << ", should be two."
            << exit(FatalError);
    }

    List<word> productMoleculesIon2
    (propsDict_.lookup("productsOfIonisedAtom"));

    if (productMoleculesIon2.size() != 2)
    {
        FatalErrorIn("moleculeAtomDissociationIonisation::setProperties()")
            << "Number of atom ionisation products is "
            << productMoleculesIon2.size()
            << ", should be two."
            << exit(FatalError);
    }


    productIdsDiss_.setSize(productMoleculesDiss.size());

    forAll(productMoleculesDiss, r)
    {
        if (productIdsDiss_.size() != 2)
        {
            FatalErrorIn("moleculeAtomDissociationIonisation::setProperties()")
                << "There should be two products (for the "
                << "dissociating molecule"
                << reactantIds_[r] << "), instead of "
                << productIdsDiss_.size() << nl
                << exit(FatalError);
        }

        forAll(productIdsDiss_, r)
        {
            productIdsDiss_[r] =
                findIndex(cloud_.typeIdList(), productMoleculesDiss[r]);

            // check that reactants belong to the typeIdList
            if (productIdsDiss_[r] == -1)
            {
                FatalErrorIn("moleculeAtomDissociationIonisation::setProperties()")
                    << "Cannot find type id: " << productMoleculesDiss[r]
                    << nl
                    << exit(FatalError);
            }
        }

        // check that product one is an 'ATOM' (not a 'MOLECULE')

        scalar rDof3 =
            cloud_.constProps(productIdsDiss_[0]).rotationalDegreesOfFreedom();

        if (rDof3 != 0)
        {
            FatalErrorIn("moleculeAtomDissociationIonisation::setProperties()")
                << "First dissociation product must be "
                << "an atom (not a molecule): " << productMoleculesDiss[0]
                << nl
                << exit(FatalError);
        }

        // check that product two is an 'ATOM' (not a 'MOLECULE')

        scalar rDof4 =
            cloud_.constProps(productIdsDiss_[1]).rotationalDegreesOfFreedom();

        if (rDof4 != 0)
        {
            FatalErrorIn("moleculeAtomDissociationIonisation::setProperties()")
                << "Second dissociation product must be "
                << "an atom (not a molecule): " << productMoleculesDiss[1]
                << nl
                << exit(FatalError);
        }
    }


    productIdsIon_.setSize(productMoleculesIon.size());

    forAll(productMoleculesIon, r)
    {
        if (productIdsIon_.size() != 2)
        {
            FatalErrorIn("moleculeAtomDissociationIonisation::setProperties()")
                << "There should be two products (for the ionising molecule "
                << reactantIds_[r] << "), instead of "
                << productIdsIon_.size() << nl
                << exit(FatalError);
        }

        forAll(productIdsIon_, r)
        {
            productIdsIon_[r] =
                findIndex(cloud_.typeIdList(), productMoleculesIon[r]);

            // check that reactants belong to the typeIdList
            if (productIdsIon_[r] == -1)
            {
                FatalErrorIn("moleculeAtomDissociationIonisation::setProperties()")
                    << "Cannot find type id: " << productMoleculesIon[r] << nl
                    << exit(FatalError);
            }
        }

        // check that product one is an 'MOLECULE' (not an 'ATOM')

        scalar rDof5 =
            cloud_.constProps(productIdsIon_[0]).rotationalDegreesOfFreedom();

        if (rDof5 < 0)
        {
            FatalErrorIn("moleculeAtomDissociationIonisation::setProperties()")
                << "First ionisation product must be a "
                << "charged molecule (not an atom/electron): "
                << productMoleculesIon[0]
                << nl
                << exit(FatalError);
        }

        // check that product two is a 'ELECTRON'

        const label charge = cloud_.constProps(productIdsIon_[1]).charge();

        if (charge != -1)
        {
            FatalErrorIn("moleculeAtomDissociationIonisation::setProperties()")
                << "Second ionisation product must be an electron: "
                << productMoleculesIon[1]
                << nl
                << exit(FatalError);
        }
    }


    productIdsIon2_.setSize(productMoleculesIon2.size());

    forAll(productMoleculesIon2, r)
    {
        if (productIdsIon2_.size() != 2)
        {
            FatalErrorIn("moleculeAtomDissociationIonisation::setProperties()")
                << "There should be two products (for the ionising molecule "
                << reactantIds_[r] << "), instead of "
                << productIdsIon2_.size() << nl
                << exit(FatalError);
        }

        forAll(productIdsIon2_, r)
        {
            productIdsIon2_[r] =
                findIndex(cloud_.typeIdList(), productMoleculesIon2[r]);

            // check that reactants belong to the typeIdList
            if (productIdsIon2_[r] == -1)
            {
                FatalErrorIn("moleculeAtomDissociationIonisation::setProperties()")
                    << "Cannot find type id: " << productMoleculesIon2[r]
                    << nl
                    << exit(FatalError);
            }
        }

        // check that product one is an 'ATOM'

        scalar rDof6 =
            cloud_.constProps(productIdsIon2_[0]).rotationalDegreesOfFreedom();

        if (rDof6 > 0)
        {
            FatalErrorIn("moleculeAtomDissociationIonisation::setProperties()")
                << "First ionisation product must be a charged atom: "
                << productMoleculesIon[0]
                << nl
                << exit(FatalError);
        }

        // check that product two is a 'ELECTRON'

        const label charge = cloud_.constProps(productIdsIon2_[1]).charge();

        if (charge != -1)
        {
            FatalErrorIn("moleculeAtomDissociationIonisation::setProperties()")
                << "Second ionisation product must be an electron: "
                << productMoleculesIon[1]
                << nl
                << exit(FatalError);
        }
    }
    
    heatOfReactionDiss_ =
        readScalar(propsDict_.lookup("heatOfReactionDissociation"));
    heatOfReactionIon_ =
        readScalar(propsDict_.lookup("heatOfReactionIonisationMolecule"));
    heatOfReactionIon2_ =
        readScalar(propsDict_.lookup("heatOfReactionIonisationAtom"));
}


bool moleculeAtomDissociationIonisation::
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


void moleculeAtomDissociationIonisation::reaction
(
    dsmcParcel& p,
    dsmcParcel& q,
    const DynamicList<label>& candidateList,
    const List<DynamicList<label> >& candidateSubList,
    const label candidateP,
    const List<label>& whichSubCell
)
{}


void moleculeAtomDissociationIonisation::reaction
(
    dsmcParcel& p,
    dsmcParcel& q
)
{
    label typeIdP = p.typeId();
    label typeIdQ = q.typeId();
    
//     reactionProperties(cloud_, p, q) = reactionProperties(cloud_, p, q);

    if (typeIdP == reactantIds_[0] && typeIdQ == reactantIds_[1])
    {
        scalar totalReactionProbability = 0.0;
        scalarList reactionProbabilities(3, 0.0);

        relax_ = true;

        bool dissocReaction = false;
        bool ionisationReaction = false;
        bool atomIonisationReaction = false;

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

        ionisationEnergy = cloud_.constProps(typeIdQ).ionisationTemperature()
                                *physicoChemical::k.value();

        // calculate if an ionisation of species Q is possible
        EcPPIon = reactionProperties(cloud_, p, q).translationalEnergy() + reactionProperties(cloud_, p, q).EEleQ();

        if ((EcPPIon - ionisationEnergy) > VSMALL)
        {
            totalReactionProbability += 1.0;
            reactionProbabilities[1] = 1.0;
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
            reactionProbabilities[2] = 1.0;
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
                            // Atom ionisation is to occur
                            atomIonisationReaction = true;
                            break;
                        }
                        if (i == 2)
                        {
                            // Dissociation reaction is to occur
                            dissocReaction = true;
                            break;
                        }
                    }
                }
            }
        }

        if (dissocReaction)
        {
            nDissociationReactionsPerTimeStep_++;
            nTotReactionsDiss_++;

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
        if (ionisationReaction)
        {
            nIonisationReactionsPerTimeStep_++;
            nTotReactionsIon_++;

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
        if (atomIonisationReaction)
        {
            nIonisationReactions2PerTimeStep_++;
            nTotReactionsIon2_++;

            if (allowSplitting_)
            {
                relax_ = false;
                
                reactionProperties(cloud_, p, q).ioniseQ
                (
                    heatOfReactionIon2_,
                    productIdsIon2_,
                    p,
                    q
                );
            }
        }
    }

    if (typeIdP == reactantIds_[1] && typeIdQ == reactantIds_[0])
    {
        scalar totalReactionProbability = 0.0;
        scalarList reactionProbabilities(3, 0.0);

        relax_ = true;

        bool dissocReaction = false;
        bool ionisationReaction = false;
        bool atomIonisationReaction = false;

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

        ionisationEnergy = cloud_.constProps(typeIdP).ionisationTemperature()
                                *physicoChemical::k.value();

        // calculate if an ionisation of species P is possible
        EcPQIon = reactionProperties(cloud_, p, q).translationalEnergy() + reactionProperties(cloud_, p, q).EEleP();

        if ((EcPQIon - ionisationEnergy) > VSMALL)
        {
            totalReactionProbability += 1.0;
            reactionProbabilities[1] = 1.0;
        }

        scalar EcPQDiss = 0.0;
        label idQ = reactionProperties(cloud_, p, q).charDissLevelQ();
        label imaxQ = 0;

        // calculate if a dissociation of species Q is possible
        EcPQDiss = reactionProperties(cloud_, p, q).translationalEnergy() + reactionProperties(cloud_, p, q).EVibQ();

        imaxQ = EcPQDiss/(physicoChemical::k.value()*reactionProperties(cloud_, p, q).thetaVQ());

        if (imaxQ - idQ > 0)
        {
           // Dissociation can occur
            totalReactionProbability += 1.0;
            reactionProbabilities[2] = 1.0;
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
                            // Atom ionisation reaction is to occur
                            atomIonisationReaction = true;
                            break;
                        }
                        if (i == 2)
                        {
                            // Dissociation reaction is to occur
                            dissocReaction = true;
                            break;
                        }
                    }
                }
            }
        }

        if (dissocReaction)
        {
            nDissociationReactionsPerTimeStep_++;
            nTotReactionsDiss_++;

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

        if (ionisationReaction)
        {
            nIonisationReactionsPerTimeStep_++;
            nTotReactionsIon_++;

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

        if (atomIonisationReaction)
        {
            nIonisationReactions2PerTimeStep_++;
            nTotReactionsIon2_++;

            if (allowSplitting_)
            {
                relax_ = false;
                
                reactionProperties(cloud_, p, q).ioniseP
                (
                    heatOfReactionIon2_,
                    productIdsIon2_,
                    p,
                    q
                );
            }
        }
    }
}


void  moleculeAtomDissociationIonisation::
outputResults(const label counterIndex)
{
    if (writeRatesToTerminal_ == true)
    {
        const List<DynamicList<dsmcParcel*> >& cellOccupancy
            = cloud_.cellOccupancy();

        List<label> mols(2, 0);
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
        label nTotReactionsDiss = nTotReactionsDiss_;
        label nTotReactionsIon = nTotReactionsIon_;
        label nTotReactionsIon2 = nTotReactionsIon2_;

        // Parallel communication
        if (Pstream::parRun())
        {
            reduce(volume, sumOp<scalar>());
            reduce(mols[0], sumOp<label>());
            reduce(mols[1], sumOp<label>());
            reduce(nTotReactionsDiss, sumOp<label>());
            reduce(nTotReactionsIon, sumOp<label>());
            reduce(nTotReactionsIon2, sumOp<label>());
        }

        numberDensities_[0] = (mols[0]*cloud().nParticle())/volume;
        numberDensities_[1] = (mols[1]*cloud().nParticle())/volume;

        const scalar deltaT = mesh_.time().deltaT().value();

        const word& reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
        const word& reactantMolB = cloud_.typeIdList()[reactantIds_[1]];

        const word& productMolA = cloud_.typeIdList()[productIdsDiss_[0]];
        const word& productMolB = cloud_.typeIdList()[productIdsDiss_[1]];

        const word& productMolC = cloud_.typeIdList()[productIdsIon_[0]];
        const word& productMolD = cloud_.typeIdList()[productIdsIon_[1]];

        const word& productMolE = cloud_.typeIdList()[productIdsIon2_[0]];
        const word& productMolF = cloud_.typeIdList()[productIdsIon2_[1]];

        if ((numberDensities_[0] > 0.0) && (numberDensities_[1] > 0.0))
        {
            scalar reactionRateDiss = 0.0;
            scalar reactionRateIon = 0.0;
            scalar reactionRateIon2 = 0.0;

            reactionRateDiss =
            (
                nTotReactionsDiss
               *cloud_.nParticle()
            )
           /(
                counterIndex*deltaT*numberDensities_[0]
               *numberDensities_[1]*volume
            );

            Info<< "Dissociation type II reaction "
                <<  reactantMolA << " + " << reactantMolB
                <<  " --> "
                << productMolA << " + " << productMolB << " + " << reactantMolB
                << ", reaction rate = " << reactionRateDiss
                << endl;

            reactionRateIon =
            (
                nTotReactionsIon
               *cloud_.nParticle()
            )
           /(
                counterIndex*deltaT*numberDensities_[0]
               *numberDensities_[1]*volume
            );

            Info<< "Ionisation reaction "
                <<  reactantMolA << " + " << reactantMolB
                <<  " --> "
                << productMolC << " + " << productMolD << " + " << reactantMolB
                << ", reaction rate = " << reactionRateIon
                << endl;

            reactionRateIon2 =
                (
                    nTotReactionsIon2
                    * cloud_.nParticle()
                )
               /(
                    counterIndex*deltaT*numberDensities_[0]
                   *numberDensities_[1]*volume
                );

            Info<< "Ionisation reaction "
                <<  reactantMolA << " + " << reactantMolB
                <<  " --> "
                << reactantMolA << " + " << productMolE << " + " << productMolF
                << ", reaction rate = " << reactionRateIon2
                << endl;
        }
    }
    else
    {
        label nTotReactionsDiss = nTotReactionsDiss_;
        label nTotReactionsIon = nTotReactionsIon_;
        label nTotReactionsIon2 = nTotReactionsIon2_;

        label nDissociationReactionsPerTimeStep =
            nDissociationReactionsPerTimeStep_;
        label nIonisationReactionsPerTimeStep =
            nIonisationReactionsPerTimeStep_;
        label nIonisationReactions2PerTimeStep =
            nIonisationReactions2PerTimeStep_;

        if (Pstream::parRun())
        {
            reduce(nTotReactionsDiss, sumOp<label>());
            reduce(nTotReactionsIon, sumOp<label>());
            reduce(nTotReactionsIon2, sumOp<label>());

            reduce(nDissociationReactionsPerTimeStep, sumOp<label>());
            reduce(nIonisationReactionsPerTimeStep, sumOp<label>());
            reduce(nIonisationReactions2PerTimeStep, sumOp<label>());
        }

        if (nTotReactionsDiss > VSMALL)
        {
            const word& reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
            const word& reactantMolB = cloud_.typeIdList()[reactantIds_[1]];

            const word& productMolA = cloud_.typeIdList()[productIdsDiss_[0]];
            const word& productMolB = cloud_.typeIdList()[productIdsDiss_[1]];

            Info<< "Dissociation type II reaction "
                <<  reactantMolA << " + " << reactantMolB
                <<  " --> "
                << productMolA << " + " << productMolB << " + "
                << reactantMolB
                << " is active, nReactions this time step = "
                << nDissociationReactionsPerTimeStep << endl;
        }

        if (nTotReactionsIon > VSMALL)
        {
            const word& reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
            const word& reactantMolB = cloud_.typeIdList()[reactantIds_[1]];

            const word& productMolA = cloud_.typeIdList()[productIdsIon_[0]];
            const word& productMolB = cloud_.typeIdList()[productIdsIon_[1]];

            Info<< "Ionisation reaction "
                <<  reactantMolA << " + " << reactantMolB
                <<  " --> "
                << productMolA << " + " << productMolB << " + "
                << reactantMolB
                << " is active, nReactions this time step = "
                << nIonisationReactionsPerTimeStep << endl;
        }

        if (nTotReactionsIon2 > VSMALL)
        {
            const word& reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
            const word& reactantMolB = cloud_.typeIdList()[reactantIds_[1]];

            const word& productMolA = cloud_.typeIdList()[productIdsIon2_[0]];
            const word& productMolB = cloud_.typeIdList()[productIdsIon2_[1]];

            Info<< "Ionisation reaction "
                <<  reactantMolA << " + " << reactantMolB
                <<  " --> "
                << reactantMolA << " + " << productMolA << " + "
                << productMolB
                << " is active, nReactions this time step = "
                << nIonisationReactions2PerTimeStep << endl;
        }
    }

    nDissociationReactionsPerTimeStep_ = 0.0;
    nIonisationReactionsPerTimeStep_ = 0.0;
    nIonisationReactions2PerTimeStep_ = 0.0;
}


bool moleculeAtomDissociationIonisation::relax() const
{
    return relax_;
}

} // End namespace Foam

// ************************************************************************* //
