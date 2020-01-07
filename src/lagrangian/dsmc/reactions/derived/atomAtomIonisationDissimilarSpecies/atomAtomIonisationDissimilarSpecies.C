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

#include "atomAtomIonisationDissimilarSpecies.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(atomAtomIonisationDissimilarSpecies, 0);

addToRunTimeSelectionTable
(
    dsmcReaction,
    atomAtomIonisationDissimilarSpecies,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

atomAtomIonisationDissimilarSpecies::atomAtomIonisationDissimilarSpecies
(
    Time& t,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcReaction(t, cloud, dict),
    productIdsIon_(),
    chargedAtom_(false),
    heatOfReactionIon_(),
    heatOfReactionIon2_(),
    nTotIonisationReactions_(0),
    nTotIonisationReactions2_(0),
    nIonisationReactionsPerTimeStep_(0),
    nIonisationReactions2PerTimeStep_(0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

atomAtomIonisationDissimilarSpecies::~atomAtomIonisationDissimilarSpecies()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void atomAtomIonisationDissimilarSpecies::initialConfiguration()
{
    setCommonReactionProperties();
    setProperties();
}


void atomAtomIonisationDissimilarSpecies::setProperties()
{

    if (reactants_[0] == reactants_[1])
    {
        FatalErrorIn("atomAtomIonisationDissimilarSpecies::setProperties()")
            << "Reactant molecules cannot be same species." << nl
            << exit(FatalError);
    }

    chargedAtom_ = Switch(propsDict_.lookup("chargedAtom"));

    // check that reactant one is an 'ATOM'

    if (rDof1_ > VSMALL)
    {
        FatalErrorIn("atomAtomIonisationDissimilarSpecies::setProperties()")
            << "First reactant must be an atom "
            << "(not a molecule or an electron): " << reactants_[0]
            << nl
            << exit(FatalError);
    }
    
    if (vDof1_ > VSMALL)
    {
         FatalErrorIn("atomAtomIonisationDissimilarSpecies::setProperties()")
            << "Reactions are currently only implemented for "
            << "monatomic and diatomic species"
            << " This is a polyatomic:" << reactants_[0]
            << nl
            << exit(FatalError);
    }

    // check that reactant two is an 'ATOM'

    if (rDof2_ > VSMALL)
    {
        FatalErrorIn("atomAtomIonisationDissimilarSpecies::setProperties()")
            << "Second reactant must be an atom "
            << "(not a molecule or an electron): " << reactants_[1]
            << nl
            << exit(FatalError);
    }

    if (vDof2_ > VSMALL)
    {
         FatalErrorIn("atomAtomIonisationDissimilarSpecies::setProperties()")
            << "Reactions are currently only implemented for "
            << "monatomic and diatomic species"
            << " This is a polyatomic:" << reactants_[1]
            << nl
            << exit(FatalError);
    }
    
    //check that reactant one is not an ion or an electron
    
    if (charge1_ == -1)
    {
        FatalErrorIn("atomAtomIonisationDissimilarSpecies::setProperties()")
            << "First reactant cannot be an ion or an electron"
            << reactants_[0]
            << nl
            << exit(FatalError);
    }

    //check that reactant two is not an ion or an electron

    if (charge2_ == -1)
    {
        FatalErrorIn("atomAtomIonisationDissimilarSpecies::setProperties()")
            << "Second reactant cannot be an ion or an electron"
            << reactants_[1]
            << nl
            << exit(FatalError);
    }

    if (chargedAtom_)
    {
        // reading in ionisation products

        List<word> productMoleculesIonisation
                            (propsDict_.lookup("productsOfIonisedAtom"));

        if (productMoleculesIonisation.size() != 2)
        {
            FatalErrorIn("atomAtomIonisationDissimilarSpecies::setProperties()")
                << "Number of ionisation products is "
                << productMoleculesIonisation.size()
                << ", should be two."
                << exit(FatalError);
        }


        productIdsIon_.setSize(productMoleculesIonisation.size());

        forAll(productMoleculesIonisation, r)
        {
            if (productIdsIon_.size() != 2)
            {
                FatalErrorIn("atomAtomIonisationDissimilarSpecies::setProperties()")
                    << "There should be two products (for the ionising "
                    << "molecule "
                    << reactants_[r] << "), instead of "
                    << productIdsIon_.size() << nl
                    << exit(FatalError);
            }

            forAll(productIdsIon_, r)
            {
                productIdsIon_[r] = findIndex(cloud_.typeIdList(),
                                              productMoleculesIonisation[r]);

                // check that reactants belong to the typeIdList
                if (productIdsIon_[r] == -1)
                {
                    FatalErrorIn("atomAtomIonisationDissimilarSpecies::setProperties()")
                        << "Cannot find type id: "
                        << productMoleculesIonisation[r] << nl
                        << exit(FatalError);
                }
            }

            // check that product three is a 'ATOM'

            scalar rDof3
            =cloud_.constProps(productIdsIon_[0]).rotationalDegreesOfFreedom();
            
            scalar vDof3
            =cloud_.constProps(productIdsIon_[0]).vibrationalDegreesOfFreedom();

            if (rDof3 > 1)
            {
                FatalErrorIn("atomAtomIonisationDissimilarSpecies::setProperties()")
                    << "First product must be an atom (not a molecule): "
                    << productMoleculesIonisation[0]
                    << nl
                    << exit(FatalError);
            }
            
            if (vDof3 > VSMALL)
            {
                FatalErrorIn("atomAtomIonisationDissimilarSpecies::setProperties()")
                    << "Reactions are currently only implemented for "
                    << "monatomic and diatomic species"
                    << " This is a polyatomic:" << productIdsIon_[1]
                    << nl
                    << exit(FatalError);
            }

            // check that product two is an 'ELECTRON'

            const label charge = cloud_.constProps(productIdsIon_[1]).charge();

            if (charge != -1)
            {
                FatalErrorIn("atomAtomIonisationDissimilarSpecies::setProperties()")
                    << "Second product must be an electron: "
                    << productMoleculesIonisation[1]
                    << nl
                    << exit(FatalError);
            }
        }
    }
    else
    {
        // reading in ionisation products

        List<word> productMoleculesIonisation
                        (propsDict_.lookup("productsOfIonisedAtomP"));

        if (productMoleculesIonisation.size() != 2)
        {
            FatalErrorIn("atomAtomIonisationDissimilarSpecies::setProperties()")
                << "Number of ionisation products is "
                << productMoleculesIonisation.size()
                << ", should be two."
                << exit(FatalError);
        }


        productIdsIon_.setSize(productMoleculesIonisation.size());

        forAll(productMoleculesIonisation, r)
        {
            if (productIdsIon_.size() != 2)
            {
                FatalErrorIn("atomAtomIonisationDissimilarSpecies::setProperties()")
                    << "There should be two products "
                    << "(for the ionising molecule "
                    << reactants_[r] << "), instead of "
                    << productIdsIon_.size() << nl
                    << exit(FatalError);
            }

            forAll(productIdsIon_, r)
            {
                productIdsIon_[r] = findIndex(cloud_.typeIdList(),
                                              productMoleculesIonisation[r]);

                // check that reactants belong to the typeIdList
                if (productIdsIon_[r] == -1)
                {
                    FatalErrorIn("atomAtomIonisationDissimilarSpecies::setProperties()")
                        << "Cannot find type id: "
                        << productMoleculesIonisation[r] << nl
                        << exit(FatalError);
                }
            }

            // check that product one is a 'ATOM'

            scalar rDof4 =
             cloud_.constProps(productIdsIon_[0]).rotationalDegreesOfFreedom();
             
            scalar vDof4 =
             cloud_.constProps(productIdsIon_[0]).vibrationalDegreesOfFreedom();

            if (rDof4 > 1)
            {
                FatalErrorIn("atomAtomIonisationDissimilarSpecies::setProperties()")
                    << "First product must be an atom (not an atom): "
                    << productMoleculesIonisation[0]
                    << nl
                    << exit(FatalError);
            }
            
            if (vDof4 > VSMALL)
            {
                FatalErrorIn("atomAtomIonisationDissimilarSpecies::setProperties()")
                    << "Reactions are currently only implemented for "
                    << "monatomic and diatomic species"
                    << " This is a polyatomic:" 
                    << productMoleculesIonisation[0]
                    << nl
                    << exit(FatalError);
            }

            // check that product two is an 'ELECTRON'

            label charge = cloud_.constProps(productIdsIon_[1]).charge();

            if (charge != -1)
            {
                FatalErrorIn("atomAtomIonisationDissimilarSpecies::setProperties()")
                    << "Second product must be an electron: "
                    << productMoleculesIonisation[1]
                    << nl
                    << exit(FatalError);
            }
        }

        List<word> productMoleculesIonisation2
                    (propsDict_.lookup("productsOfIonisedAtomQ"));

        if (productMoleculesIonisation2.size() != 2)
        {
            FatalErrorIn("atomAtomIonisationDissimilarSpecies::setProperties()")
                << "Number of ionisation products is "
                << productMoleculesIonisation2.size()
                << ", should be two."
                << exit(FatalError);
        }


        productIdsIon2_.setSize(productMoleculesIonisation2.size());

        forAll(productMoleculesIonisation2, r)
        {
            if (productIdsIon2_.size() != 2)
            {
                FatalErrorIn("atomAtomIonisationDissimilarSpecies::setProperties()")
                    << "There should be two products "
                    << "(for the ionising molecule "
                    << reactants_[r] << "), instead of "
                    << productIdsIon2_.size() << nl
                    << exit(FatalError);
            }

            forAll(productIdsIon2_, r)
            {
                productIdsIon2_[r] = findIndex(cloud_.typeIdList(),
                                               productMoleculesIonisation2[r]);

                // check that reactants belong to the typeIdList
                if (productIdsIon2_[r] == -1)
                {
                    FatalErrorIn("atomAtomIonisationDissimilarSpecies::setProperties()")
                        << "Cannot find type id: "
                        << productMoleculesIonisation2[r] << nl
                        << exit(FatalError);
                }
            }

            // check that product one is a 'ATOM', not an 'MOLECULE'

            scalar rDof5 =
            cloud_.constProps(productIdsIon2_[0]).rotationalDegreesOfFreedom();
            
            scalar vDof5 =
            cloud_.constProps(productIdsIon2_[0]).vibrationalDegreesOfFreedom();

            if (rDof5 > 1)
            {
                FatalErrorIn("atomAtomIonisationDissimilarSpecies::setProperties()")
                    << "First product must be an atom (not an atom): "
                    << productMoleculesIonisation2[0]
                    << nl
                    << exit(FatalError);
            }
            
            if (vDof5 > VSMALL)
            {
                FatalErrorIn("atomAtomIonisationDissimilarSpecies::setProperties()")
                    << "Reactions are currently only implemented for "
                    << "monatomic and diatomic species"
                    << " This is a polyatomic:" 
                    << productMoleculesIonisation2[0]
                    << nl
                    << exit(FatalError);
            }

            // check that product two is an 'ELECTRON'

            const label charge =
                            cloud_.constProps(productIdsIon2_[1]).charge();

            if (charge != -1)
            {
                FatalErrorIn("atomAtomIonisationDissimilarSpecies::setProperties()")
                    << "Second product must be an electron: "
                    << productMoleculesIonisation2[1]
                    << nl
                    << exit(FatalError);
            }
        }
    }

    if (chargedAtom_)
    {
        heatOfReactionIon_ =
                readScalar(propsDict_.lookup("heatOfReactionIonisation"));
    }
    else
    {
        heatOfReactionIon_ =
                readScalar(propsDict_.lookup("heatOfReactionIonisationP"));
        heatOfReactionIon2_ =
            readScalar(propsDict_.lookup("heatOfReactionIonisationQ"));
    }
}


bool atomAtomIonisationDissimilarSpecies::tryReactMolecules(
                                label typeIdP, const label typeIdQ)
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


void atomAtomIonisationDissimilarSpecies::reaction
(
    dsmcParcel& p,
    dsmcParcel& q,
    const DynamicList<label>& candidateList,
    const List<DynamicList<label> >& candidateSubList,
    const label candidateP,
    const List<label>& whichSubCell
)
{}


void atomAtomIonisationDissimilarSpecies::reaction
(
    dsmcParcel& p,
    dsmcParcel& q
)
{
    label typeIdP = p.typeId();
    label typeIdQ = q.typeId();
    
//     auto rp = reactionProperties(cloud_, p, q);

    if (typeIdP == reactantIds_[0] && typeIdQ == reactantIds_[1]
        && !chargedAtom_)
    {
        relax_ = true;

        scalar totalReactionProbability = 0.0;
        scalarList reactionProbabilities(2, 0.0);

        bool ionisationReactionP = false;
        bool ionisationReactionQ = false;

        scalar EcPPIon = 0.0;
        scalar ionisationEnergy =
                        cloud_.constProps(typeIdP).ionisationTemperature()
                        *physicoChemical::k.value();

        // calculate if an ionisation of particle P is possible
        EcPPIon = reactionProperties(cloud_, p, q).translationalEnergy() + reactionProperties(cloud_, p, q).EEleP();

        if ((EcPPIon - ionisationEnergy) > VSMALL)
        {
            totalReactionProbability += 1.0;
            reactionProbabilities[0] = 1.0;
        }

        ionisationEnergy = cloud_.constProps(typeIdQ).ionisationTemperature()
                                *physicoChemical::k.value();

        // calculate if an ionisation of particle Q is possible
        EcPPIon = reactionProperties(cloud_, p, q).translationalEnergy() + reactionProperties(cloud_, p, q).EEleQ();

        if ((EcPPIon - ionisationEnergy) > VSMALL)
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
                            ionisationReactionP = true;
                            break;
                        }
                        if (i == 1)
                        {
                            ionisationReactionQ = true;
                            break;
                        }
                    }
                }
            }
        }

        if (ionisationReactionP)
        {
            nTotIonisationReactions_++;
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

        if (ionisationReactionQ)
        {
            nTotIonisationReactions2_++;
            nIonisationReactions2PerTimeStep_++;

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

    if (typeIdP == reactantIds_[1] && typeIdQ == reactantIds_[0]
        && !chargedAtom_)
    {
        relax_ = true;

        scalar totalReactionProbability = 0.0;
        scalarList reactionProbabilities(2, 0.0);

        bool ionisationReactionP = false;
        bool ionisationReactionQ = false;

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
                            ionisationReactionQ = true;
                            break;
                        }
                        if (i == 1)
                        {
                            ionisationReactionP = true;
                            break;
                        }
                    }
                }
            }
        }

        if (ionisationReactionP)
        {
            nTotIonisationReactions_++;
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

        if (ionisationReactionQ)
        {
            nTotIonisationReactions2_++;
            nIonisationReactions2PerTimeStep_++;

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

    if (typeIdP == reactantIds_[0] && typeIdQ == reactantIds_[1]
        && chargedAtom_)
    {
        relax_ = true;

        scalar ionisationEnergy =
                    cloud_.constProps(typeIdP).ionisationTemperature()
                    *physicoChemical::k.value();

        // calculate if an ionisation of species P is possible
        scalar EcPPIon = reactionProperties(cloud_, p, q).translationalEnergy() + reactionProperties(cloud_, p, q).EEleP();

        if ((EcPPIon - ionisationEnergy) > VSMALL)
        {
            // P can ionise
            nTotIonisationReactions_++;
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

    if (typeIdP == reactantIds_[1] && typeIdQ == reactantIds_[0]
        && chargedAtom_)
    {
        relax_ = true;

        // calculate if an ionisation of species Q is possible
        scalar EcPPIon = reactionProperties(cloud_, p, q).translationalEnergy() + reactionProperties(cloud_, p, q).EEleQ();

        scalar ionisationEnergy =
            cloud_.constProps(typeIdQ).ionisationTemperature()
           *physicoChemical::k.value();

        if ((EcPPIon - ionisationEnergy) > VSMALL)
        {
            // Q can ionise

            nTotIonisationReactions_++;
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


void  atomAtomIonisationDissimilarSpecies::outputResults(const label
                                                    counterIndex)
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
        label nTotReactionsIonisation = nTotIonisationReactions_;
        label nTotReactionsIonisation2 = nTotIonisationReactions2_;

        // Parallel communication
        if (Pstream::parRun())
        {
            reduce(volume, sumOp<scalar>());
            reduce(mols[0], sumOp<label>());
            reduce(mols[1], sumOp<label>());
            reduce(nTotReactionsIonisation, sumOp<label>());
            reduce(nTotReactionsIonisation2, sumOp<label>());
        }

        numberDensities_[0] = (mols[0]*cloud().nParticle())/volume;
        numberDensities_[1] = (mols[1]*cloud().nParticle())/volume;

        const scalar deltaT = mesh_.time().deltaT().value();

        word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
        word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];

        word productMolC = cloud_.typeIdList()[productIdsIon_[0]];
        word productMolD = cloud_.typeIdList()[productIdsIon_[1]];

        word productMolE;
        word productMolF;

        if (!chargedAtom_)
        {
            productMolE = cloud_.typeIdList()[productIdsIon2_[0]];
            productMolF = cloud_.typeIdList()[productIdsIon2_[1]];
        }

        if ((numberDensities_[0] > 0.0) && (numberDensities_[1] > 0.0))
        {
            scalar reactionRateIonisation = 0.0;
            scalar reactionRateIonisation2 = 0.0;

            reactionRateIonisation =
            (
            nTotReactionsIonisation
            * cloud_.nParticle()
)/(counterIndex*deltaT*numberDensities_[0]
                    * numberDensities_[1]*volume);

            Info<< "Ionisation reaction "
                <<  reactantMolA << " + " << reactantMolB
                <<  " --> "
                << productMolC << " + " << productMolD << " + " << reactantMolB
                << ", reaction rate = " << reactionRateIonisation
                << endl;

            if (!chargedAtom_)
            {
                reactionRateIonisation2 =
                (
                nTotReactionsIonisation2
                * cloud_.nParticle()
)/(counterIndex*deltaT*numberDensities_[0]
                        * numberDensities_[1]*volume);

                Info<< "Ionisation reaction "
                    <<  reactantMolA << " + " << reactantMolB
                    <<  " --> "
                    << reactantMolA << " + " << productMolE
                    << " + " << productMolF
                    << ", reaction rate = " << reactionRateIonisation2
                    << endl;
            }
        }
    }
    else
    {
        label nTotReactionsIonisation = nTotIonisationReactions_;
        label nTotReactionsIonisation2 = nTotIonisationReactions2_;
        label nIonisationReactionsPerTimeStep =
                    nIonisationReactionsPerTimeStep_;
        label nIonisationReactions2PerTimeStep =
                    nIonisationReactions2PerTimeStep_;

        if (Pstream::parRun())
        {
            reduce(nTotReactionsIonisation, sumOp<label>());
            reduce(nTotReactionsIonisation2, sumOp<label>());
            reduce(nIonisationReactionsPerTimeStep, sumOp<label>());
            reduce(nIonisationReactions2PerTimeStep, sumOp<label>());
        }

        if (nTotReactionsIonisation > VSMALL)
        {
                word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
                word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];

                word productMolA = cloud_.typeIdList()[productIdsIon_[0]];
                word productMolB = cloud_.typeIdList()[productIdsIon_[1]];

                Info<< "Ionisation reaction "
                <<  reactantMolA << " + " << reactantMolB
                <<  " --> "
                << productMolA << " + " << productMolB
                << " + " << reactantMolB
                << " is active, nReactions this time step = "
                << nIonisationReactionsPerTimeStep << endl;
        }

        if (nTotReactionsIonisation2 > VSMALL)
        {
                word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
                word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];

                word productMolA = cloud_.typeIdList()[productIdsIon2_[0]];
                word productMolB = cloud_.typeIdList()[productIdsIon2_[1]];

                Info<< "Ionisation reaction "
                <<  reactantMolA << " + " << reactantMolB
                <<  " --> "
                << reactantMolA << " + " << productMolA
                << " + " << productMolB
                << " is active, nReactions this time step = "
                << nIonisationReactions2PerTimeStep << endl;
        }
    }

    nIonisationReactionsPerTimeStep_ = 0.0;
    nIonisationReactions2PerTimeStep_ = 0.0;
}


bool atomAtomIonisationDissimilarSpecies::relax() const
{
    return relax_;
}

} // End namespace Foam

// ************************************************************************* //
