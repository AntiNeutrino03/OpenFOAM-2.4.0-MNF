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

#include "forwardAssociativeIonisationDissimilarSpecies.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(forwardAssociativeIonisationDissimilarSpecies, 0);

addToRunTimeSelectionTable
(
    dsmcReaction,
    forwardAssociativeIonisationDissimilarSpecies,
    dictionary
);




// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
forwardAssociativeIonisationDissimilarSpecies::
    forwardAssociativeIonisationDissimilarSpecies
(
    Time& t,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcReaction(t, cloud, dict),
    intermediateId_(),
    associativeIonisationProductIds_(),
    ionisationPProductIds_(),
    ionisationQProductIds_(),
    heatOfReactionRecombination_(),
    heatOfReactionIntermediateIonisation_(),
    heatOfReactionIonisationP_(),
    heatOfReactionIonisationQ_(),
    nTotalIonisationPReactions_(0),
    nTotalIonisationQReactions_(0),
    nTotalAssociativeIonisationReactions_(0),
    nIonisationPReactionsPerTimeStep_(0),
    nIonisationQReactionsPerTimeStep_(0),
    nAssociativeIonisationReactionsPerTimeStep_(0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

forwardAssociativeIonisationDissimilarSpecies::
~forwardAssociativeIonisationDissimilarSpecies()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void forwardAssociativeIonisationDissimilarSpecies::initialConfiguration()
{
    setCommonReactionProperties();
    setProperties();
}


void forwardAssociativeIonisationDissimilarSpecies::setProperties()
{
    forAll(reactantIds_, r)
    {
        // check that reactants are 'ATOMS' (not 'MOLECULES')

        scalar rDofReactant =
            cloud_.constProps(reactantIds_[r]).rotationalDegreesOfFreedom();

        if (rDofReactant > VSMALL)
        {
            FatalErrorIn("forwardAssociativeIonisationDissimilarSpecies::setProperties()")
                << "Reactant must be an atom (not a molecule): "
                << reactants_[r]
                << nl
                << exit(FatalError);
        }
    }

    // reading in associative ionisation products

    const List<word> associativeIonisationProductMolecules
    (propsDict_.lookup("productsOfAssociativeIonisation"));

    associativeIonisationProductIds_.setSize
    (associativeIonisationProductMolecules.size(), -1);

    forAll(associativeIonisationProductIds_, i)
    {
        associativeIonisationProductIds_[i] = findIndex(cloud_.typeIdList(),
                                    associativeIonisationProductMolecules[i]);

        // check that products belong to the typeIdList
        if (associativeIonisationProductIds_[i] == -1)
        {
            FatalErrorIn("forwardAssociativeIonisationDissimilarSpecies::setProperties()")
                << "Cannot find type id: "
                << associativeIonisationProductMolecules[i] << nl
                << exit(FatalError);
        }
    }


    // check that products are a 'MOLECULE' and an 'ELECTRON'

    label id0 = associativeIonisationProductIds_[0];

    scalar rDofProd1 = cloud_.constProps(id0).rotationalDegreesOfFreedom();

    if (rDofProd1 < 1)
    {
        FatalErrorIn("forwardAssociativeIonisationDissimilarSpecies::setProperties()")
            << "First product must be a molecule: "
            << associativeIonisationProductMolecules
            << nl
            << exit(FatalError);
    }

    label id1 = associativeIonisationProductIds_[1];

    label charge = cloud_.constProps(id1).charge();

    if (charge != -1)
    {
        FatalErrorIn("forwardAssociativeIonisationDissimilarSpecies::setProperties()")
            << "Second product must be an electron: "
            << associativeIonisationProductMolecules
            << nl
            << exit(FatalError);
    }

    // reading in ionisation of P products

    const List<word> ionisationPProductMolecules
    (propsDict_.lookup("productsOfIonisationP"));

    ionisationPProductIds_.setSize(ionisationPProductMolecules.size(), -1);

    forAll(ionisationPProductIds_, i)
    {
        ionisationPProductIds_[i] = findIndex(cloud_.typeIdList(),
                                              ionisationPProductMolecules[i]);

        // check that products belong to the typeIdList
        if (ionisationPProductIds_[i] == -1)
        {
            FatalErrorIn("forwardAssociativeIonisationDissimilarSpecies::setProperties()")
                << "Cannot find type id: "
                << ionisationPProductMolecules[i] << nl
                << exit(FatalError);
        }
    }


    // check that products are an 'ATOM' and an 'ELECTRON'

    rDofProd1 =
     cloud_.constProps(ionisationPProductIds_[0]).rotationalDegreesOfFreedom();

    if (rDofProd1 > 0)
    {
        FatalErrorIn("forwardAssociativeIonisationDissimilarSpecies::setProperties()")
            << "First product must be a charged atom: "
            << ionisationPProductMolecules
            << nl
            << exit(FatalError);
    }

    const label charge2 =
        cloud_.constProps(ionisationPProductIds_[1]).charge();

    if (charge2 != -1)
    {
        FatalErrorIn("forwardAssociativeIonisationDissimilarSpecies::setProperties()")
            << "Second product must be an electron: "
            << ionisationPProductMolecules
            << nl
            << exit(FatalError);
    }

    // reading in ionisation of Q products

    const List<word> ionisationQProductMolecules
    (propsDict_.lookup("productsOfIonisationQ"));

    ionisationQProductIds_.setSize(ionisationQProductMolecules.size(), -1);

    forAll(ionisationQProductIds_, i)
    {
        ionisationQProductIds_[i] = findIndex(cloud_.typeIdList(),
                                              ionisationQProductMolecules[i]);

        // check that products belong to the typeIdList
        if (ionisationQProductIds_[i] == -1)
        {
            FatalErrorIn("forwardAssociativeIonisationDissimilarSpecies::setProperties()")
                << "Cannot find type id: "
                << ionisationQProductMolecules[i] << nl
                << exit(FatalError);
        }
    }


    // check that products are an 'ATOM' and an 'ELECTRON'

    rDofProd1 =
     cloud_.constProps(ionisationQProductIds_[0]).rotationalDegreesOfFreedom();

    if (rDofProd1 > 0)
    {
        FatalErrorIn("forwardAssociativeIonisationDissimilarSpecies::setProperties()")
            << "First product must be a charged atom: "
            << ionisationQProductMolecules
            << nl
            << exit(FatalError);
    }

    const label charge3 =
            cloud_.constProps(ionisationQProductIds_[1]).charge();

    if (charge3 != -1)
    {
        FatalErrorIn("forwardAssociativeIonisationDissimilarSpecies::setProperties()")
            << "Second product must be an electron: "
            << ionisationQProductMolecules
            << nl
            << exit(FatalError);
    }

    // reading in intermediate molecule

    const word intermediateMolecule
            (propsDict_.lookup("intermediateMolecule"));

    intermediateId_ = findIndex(cloud_.typeIdList(), intermediateMolecule);

    // check that reactants belong to the typeIdList (constant/dsmcProperties)
    if (intermediateId_ == -1)
    {
        FatalErrorIn("forwardAssociativeIonisationDissimilarSpecies::setProperties()")
            << "Cannot find type id: "
            << intermediateMolecule << nl
            << exit(FatalError);
    }

    // check that the intermediate is a 'MOLECULE'

    const scalar rDof =
            cloud_.constProps(intermediateId_).rotationalDegreesOfFreedom();

    if (rDof < 1)
    {
        FatalErrorIn("forwardAssociativeIonisationDissimilarSpecies::setProperties()")
            << "The intermediate specie must be a molecule (not an atom): "
            << intermediateMolecule
            << nl
            << exit(FatalError);
    }
    
    heatOfReactionRecombination_ =
        readScalar(propsDict_.lookup("heatOfReactionRecombination"));
    heatOfReactionIntermediateIonisation_ =
        readScalar(propsDict_.lookup("heatOfReactionIntermediateIonisation"));
    heatOfReactionIonisationP_ =
        readScalar(propsDict_.lookup("heatOfReactionIonisationP"));
    heatOfReactionIonisationQ_ =
        readScalar(propsDict_.lookup("heatOfReactionIonisationQ"));
}


bool forwardAssociativeIonisationDissimilarSpecies::
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


void forwardAssociativeIonisationDissimilarSpecies::reaction
(
    dsmcParcel& p,
    dsmcParcel& q,
    const DynamicList<label>& candidateList,
    const List<DynamicList<label> >& candidateSubList,
    const label candidateP,
    const List<label>& whichSubCell
)
{}


void forwardAssociativeIonisationDissimilarSpecies::reaction
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
        scalarList reactionProbabilities(3, 0.0);
        scalar omegaIntermediate =
            cloud_.constProps(intermediateId_).omega();
        scalar rotationalDofIntermediate =
            cloud_.constProps(intermediateId_).rotationalDegreesOfFreedom();
        scalar ChiBIntermediate = 2.5 - omegaIntermediate;
        scalar thetaVIntermediate =
            cloud_.constProps(intermediateId_).thetaV()[0];
        scalar thetaDIntermediate =
            cloud_.constProps(intermediateId_).thetaD()[0];
        scalar ZrefIntermediate =
            cloud_.constProps(intermediateId_).Zref()[0];
        scalar refTempZvIntermediate =
            cloud_.constProps(intermediateId_).TrefZv()[0];
        label ELevelIntermediate = 0;
        const List<scalar>& EElistIntermediate =
            cloud_.constProps(intermediateId_).electronicEnergyList();
        const List<label>& gListIntermediate =
            cloud_.constProps(intermediateId_).degeneracyList();
        label jMaxIntermediate =
            cloud_.constProps(intermediateId_).numberOfElectronicLevels();


        bool ionisationReactionP = false;
        bool ionisationReactionQ = false;
        bool associativeIonisation = false;

        // 3 reactions possible
        
        // 1. Ionisation of P
        // 2. Ionisation of Q
        // 3. Forward associative ionisation

        scalar Ec = 0.0;

        scalar ionisationEnergy =
                cloud_.constProps(typeIdP).ionisationTemperature()
                *physicoChemical::k.value();

        // calculate if an ionisation of species P is possible
        Ec = reactionProperties(cloud_, p, q).translationalEnergy() + reactionProperties(cloud_, p, q).EEleP();

        if ((Ec - ionisationEnergy) > VSMALL)
        {
            // Ionisation can occur
            totalReactionProbability += 1.0;
            reactionProbabilities[0] = 1.0;
        }

        ionisationEnergy = cloud_.constProps(typeIdQ).ionisationTemperature()
                                *physicoChemical::k.value();

        // calculate if an ionisation of species P is possible
        Ec = reactionProperties(cloud_, p, q).translationalEnergy() + reactionProperties(cloud_, p, q).EEleQ();

        if ((Ec - ionisationEnergy) > VSMALL)
        {
            // Ionisation can occur
            totalReactionProbability += 1.0;
            reactionProbabilities[1] = 1.0;
        }

        Ec = reactionProperties(cloud_, p, q).translationalEnergy() + reactionProperties(cloud_, p, q).EEleP();
        scalar EcOrig = Ec;

        label iMax = (Ec /(physicoChemical::k.value()*thetaVIntermediate));

        label vibLevelIntermediate = -1;

        if (iMax > SMALL)
        {
            vibLevelIntermediate = cloud_.postCollisionVibrationalEnergyLevel
                (
                    false,
                    0.0,
                    iMax,
                    thetaVIntermediate,
                    thetaDIntermediate,
                    refTempZvIntermediate,
                    omegaIntermediate,
                    ZrefIntermediate,
                    Ec
                    );
        }

        if (vibLevelIntermediate == 0)
        {
            //'Form' the intermediate molecule and test it for ionisation
            const scalar heatOfReactionRecombinationJoules =
                    heatOfReactionRecombination_*physicoChemical::k.value();

            Ec = EcOrig + heatOfReactionRecombinationJoules;

            iMax = (Ec /(physicoChemical::k.value()*thetaVIntermediate));

            if (iMax > SMALL)
            {
                label postCollisionVibLevel =
                cloud_.postCollisionVibrationalEnergyLevel
                (
                    true,
                    0.0,
                    iMax,
                    thetaVIntermediate,
                    thetaDIntermediate,
                    refTempZvIntermediate,
                    omegaIntermediate,
                    ZrefIntermediate,
                    Ec
                );

                Ec -= postCollisionVibLevel*thetaVIntermediate
                                *physicoChemical::k.value();
            }

            label postCollisionELevel =
                    cloud_.postCollisionElectronicEnergyLevel
                    (
                        Ec,
                        jMaxIntermediate,
                        omegaIntermediate,
                        EElistIntermediate,
                        gListIntermediate
                    );

            ELevelIntermediate = postCollisionELevel;

            // relative translational energy after electronic exchange
            Ec -= EElistIntermediate[ELevelIntermediate];

            scalar ERot = 0.0;

            scalar energyRatio = cloud_.postCollisionRotationalEnergy
                                (
                                    rotationalDofIntermediate,
                                    ChiBIntermediate
                                );


            ERot = energyRatio*Ec;

            Ec -= ERot;


            // redistribution finished, test it for ionisation

            scalar EcPPIon = 0.0;
            scalar ionisationEnergy =
                cloud_.constProps(intermediateId_).ionisationTemperature()
                *physicoChemical::k.value();

            // calculate if an ionisation of species P is possible
            EcPPIon = Ec + EElistIntermediate[ELevelIntermediate];

            if ((EcPPIon - ionisationEnergy) > VSMALL)
            {
                // Associative ionisation can occur
                totalReactionProbability += 1.0;
                reactionProbabilities[2] = 1.0;
            }
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
                            ionisationReactionP = true;
                            break;
                        }
                        if (i == 1)
                        {
                            // Ionisation reaction is to occur
                            ionisationReactionQ = true;
                            break;
                        }
                        if (i == 2)
                        {
                            // Associative ionisation is to occur
                            associativeIonisation = true;
                            break;
                        }
                    }
                }
            }
        }

        if (ionisationReactionP)
        {
            nTotalIonisationPReactions_++;
            nIonisationPReactionsPerTimeStep_++;

            if (allowSplitting_)
            {
                relax_ = false;

                reactionProperties(cloud_, p, q).ioniseP
                (
                    heatOfReactionIonisationP_,
                    ionisationPProductIds_,
                    p,
                    q
                );
            }
        }

        if (ionisationReactionQ)
        {
            nTotalIonisationQReactions_++;
            nIonisationQReactionsPerTimeStep_++;

            if (allowSplitting_)
            {
                relax_ = false;
                
                reactionProperties(cloud_, p, q).ioniseQ
                (
                    heatOfReactionIonisationQ_,
                    ionisationQProductIds_,
                    p,
                    q
                );
            }
        }

        if (associativeIonisation)
        {
            nTotalAssociativeIonisationReactions_++;
            nAssociativeIonisationReactionsPerTimeStep_++;

            if (allowSplitting_)
            {
                relax_ = false;
                
                reactionProperties(cloud_, p, q).associativeIonisation
                (
                    heatOfReactionIntermediateIonisation_,
                    heatOfReactionRecombination_,
                    associativeIonisationProductIds_,
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
        scalarList reactionProbabilities(3, 0.0);

        scalar omegaIntermediate = cloud_.constProps(intermediateId_).omega();
        scalar rotationalDofIntermediate =
            cloud_.constProps(intermediateId_).rotationalDegreesOfFreedom();
        scalar ChiBIntermediate = 2.5 - omegaIntermediate;
        scalar thetaVIntermediate =
            cloud_.constProps(intermediateId_).thetaV()[0];
        scalar thetaDIntermediate =
            cloud_.constProps(intermediateId_).thetaD()[0];
        scalar ZrefIntermediate = cloud_.constProps(intermediateId_).Zref()[0];
        scalar refTempZvIntermediate =
            cloud_.constProps(intermediateId_).TrefZv()[0];
        label ELevelIntermediate = 0;
        const List<scalar>& EElistIntermediate =
            cloud_.constProps(intermediateId_).electronicEnergyList();
        const List<label>& gListIntermediate =
            cloud_.constProps(intermediateId_).degeneracyList();
        label jMaxIntermediate =
            cloud_.constProps(intermediateId_).numberOfElectronicLevels();

        bool ionisationReactionP = false;
        bool ionisationReactionQ = false;
        bool associativeIonisation = false;

        // 3 reactions possible
        
        // 1. Ionisation of P
        // 2. Ionisation of Q
        // 3. Forward associative ionisation

        // collision energy is the translational energy of the two atoms,
        // plus their electronic energies

        scalar Ec = 0.0;

        scalar ionisationEnergy =
                cloud_.constProps(typeIdP).ionisationTemperature()
                *physicoChemical::k.value();

        // calculate if an ionisation of species Q is possible
        Ec = reactionProperties(cloud_, p, q).translationalEnergy() + reactionProperties(cloud_, p, q).EEleP();

        if ((Ec - ionisationEnergy) > VSMALL)
        {
            // Ionisation can occur
            totalReactionProbability += 1.0;
            reactionProbabilities[1] = 1.0;
        }

        ionisationEnergy = cloud_.constProps(typeIdQ).ionisationTemperature()
                            *physicoChemical::k.value();

        // calculate if an ionisation of species P is possible
        Ec = reactionProperties(cloud_, p, q).translationalEnergy() + reactionProperties(cloud_, p, q).EEleQ();

        if ((Ec - ionisationEnergy) > VSMALL)
        {
            // Ionisation can occur
            totalReactionProbability += 1.0;
            reactionProbabilities[0] = 1.0;
        }

        Ec = reactionProperties(cloud_, p, q).translationalEnergy() + reactionProperties(cloud_, p, q).EEleQ();
        scalar EcOrig = Ec;

        label iMax = (Ec /(physicoChemical::k.value()*thetaVIntermediate));
        label vibLevelIntermediate = -1;

        if (iMax > SMALL)
        {
            vibLevelIntermediate = cloud_.postCollisionVibrationalEnergyLevel
                    (
                        true,
                        0.0,
                        iMax,
                        thetaVIntermediate,
                        thetaDIntermediate,
                        refTempZvIntermediate,
                        omegaIntermediate,
                        ZrefIntermediate,
                        Ec
                    );
        }

        if (vibLevelIntermediate == 0)
        {
            //'Form' the intermediate molecule and test it for ionisation
            const scalar heatOfReactionRecombinationJoules =
                        heatOfReactionRecombination_
                        *physicoChemical::k.value();

            Ec = EcOrig + heatOfReactionRecombinationJoules;


            iMax = (Ec /(physicoChemical::k.value()*thetaVIntermediate));

            if (iMax > SMALL)
            {
                label postCollisionVibLevel =
                    cloud_.postCollisionVibrationalEnergyLevel
                    (
                        true,
                        0.0,
                        iMax,
                        thetaVIntermediate,
                        thetaDIntermediate,
                        refTempZvIntermediate,
                        omegaIntermediate,
                        ZrefIntermediate,
                        Ec
                    );

                Ec -= postCollisionVibLevel*thetaVIntermediate
                            *physicoChemical::k.value();
            }

            label postCollisionELevel =
                    cloud_.postCollisionElectronicEnergyLevel
                    (
                        Ec,
                        jMaxIntermediate,
                        omegaIntermediate,
                        EElistIntermediate,
                        gListIntermediate
                    );

            ELevelIntermediate = postCollisionELevel;

            // relative transitional energy after electronic exchange
            Ec -= EElistIntermediate[ELevelIntermediate];

            scalar ERot = 0.0;

            scalar energyRatio = cloud_.postCollisionRotationalEnergy
                                (rotationalDofIntermediate, ChiBIntermediate);

            ERot = energyRatio*Ec;

            Ec -= ERot;

            // redistribution finished, test it for ionisation

            scalar EcPPIon = 0.0;
            scalar ionisationEnergy =
                    cloud_.constProps(intermediateId_).ionisationTemperature()
                    *physicoChemical::k.value();

            // calculate if an ionisation of species P is possible
            EcPPIon = Ec + EElistIntermediate[ELevelIntermediate];

            if ((EcPPIon - ionisationEnergy) > VSMALL)
            {
                // Associative ionisation can occur
                totalReactionProbability += 1.0;
                reactionProbabilities[2] = 1.0;
            }
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
                            ionisationReactionP = true;
                            break;
                        }
                        if (i == 1)
                        {
                            // Ionisation reaction is to occur
                            ionisationReactionQ = true;
                            break;
                        }
                        if (i == 2)
                        {
                            // Associative ionisation is to occur
                            associativeIonisation = true;
                            break;
                        }
                    }
                }
            }
        }

        // Q ionises, is called P here for consistency with the reaction rates
        if (ionisationReactionP)
        {
            nTotalIonisationPReactions_++;
            nIonisationPReactionsPerTimeStep_++;

            if (allowSplitting_)
            {
                relax_ = false;
                
                reactionProperties(cloud_, p, q).ioniseQ
                (
                    heatOfReactionIonisationP_,
                    ionisationPProductIds_,
                    p,
                    q
                );
            }
        }

        // P ionises, is called Q for consistency with rates
        if (ionisationReactionQ)
        {
            nTotalIonisationQReactions_++;
            nIonisationQReactionsPerTimeStep_++;

            if (allowSplitting_)
            {
                relax_ = false;
                
                reactionProperties(cloud_, p, q).ioniseP
                (
                    heatOfReactionIonisationQ_,
                    ionisationQProductIds_,
                    p,
                    q
                );
            }
        }

        if (associativeIonisation)
        {
            nTotalAssociativeIonisationReactions_++;
            nAssociativeIonisationReactionsPerTimeStep_++;

            if (allowSplitting_)
            {
                relax_ = false;
                
                reactionProperties(cloud_, p, q).associativeIonisation
                (
                    heatOfReactionIntermediateIonisation_,
                    heatOfReactionRecombination_,
                    associativeIonisationProductIds_,
                    p,
                    q
                );
            }
        }
    }
}


void  forwardAssociativeIonisationDissimilarSpecies::
outputResults(const label counterIndex)
{
    if (writeRatesToTerminal_ == true)
    {
        // measure density

        const List<DynamicList<dsmcParcel*> >& cellOccupancy
            = cloud_.cellOccupancy();

        volume_ = 0.0;

        forAll(cellOccupancy, c)
        {
            volume_ += mesh_.cellVolumes()[c];
        }

        List<label> mols (2, 0);
        scalar volume = volume_;
        label nTotalAssociativeIonisationReactions =
                        nTotalAssociativeIonisationReactions_;
        label nTotalIonisationPReactions = nTotalIonisationPReactions_;
        label nTotalIonisationQReactions = nTotalIonisationQReactions_;

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
        }

        // Parallel communication
        if (Pstream::parRun())
        {
            reduce(mols[0], sumOp<label>());
            reduce(mols[1], sumOp<label>());
            reduce(volume, sumOp<scalar>());
            reduce(nTotalAssociativeIonisationReactions, sumOp<label>());
            reduce(nTotalIonisationPReactions, sumOp<label>());
            reduce(nTotalIonisationQReactions, sumOp<label>());
        }

        numberDensities_[0] = (mols[0]*cloud().nParticle())/volume;

        if (reactantIds_[0] == reactantIds_[1])
        {
            numberDensities_[1] = (mols[0]*cloud().nParticle())/volume;
        }
        else
        {
            numberDensities_[1] = (mols[1]*cloud().nParticle())/volume;
        }

        const scalar deltaT = mesh_.time().deltaT().value();

        word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
        word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];

        word productMolA =
            cloud_.typeIdList()[associativeIonisationProductIds_[0]];
        word productMolB =
            cloud_.typeIdList()[associativeIonisationProductIds_[1]];

        word productMolC = cloud_.typeIdList()[ionisationPProductIds_[0]];
        word productMolD = cloud_.typeIdList()[ionisationPProductIds_[1]];

        word productMolE = cloud_.typeIdList()[ionisationQProductIds_[0]];
        word productMolF = cloud_.typeIdList()[ionisationQProductIds_[1]];

        if ((numberDensities_[0] > 0.0) && (numberDensities_[1] > 0.0))
        {
            scalar reactionRateAssociativeIonisation = 0.0;
            scalar reactionRateIonisationP = 0.0;
            scalar reactionRateIonisationQ = 0.0;

            reactionRateAssociativeIonisation =
            (
                nTotalAssociativeIonisationReactions
                * cloud_.nParticle()
)/(counterIndex*deltaT*numberDensities_[0]
            * numberDensities_[1]*volume);

            Info<< "Associative ionisation reaction "
                <<  reactantMolA << " + " << reactantMolB
                <<  " --> "
                << productMolA << " + " << productMolB
                << ", reaction rate = " << reactionRateAssociativeIonisation
                << endl;

            reactionRateIonisationP =
            (
                nTotalIonisationPReactions
                * cloud_.nParticle()
)/(counterIndex*deltaT*numberDensities_[0]
            * numberDensities_[1]*volume);

            Info<< "Ionisation reaction "
                <<  reactantMolA << " + " << reactantMolB
                <<  " --> "
                << productMolC << " + " << productMolD << " + " << reactantMolB
                << ", reaction rate = " << reactionRateIonisationP
                << endl;

            reactionRateIonisationQ =
            (
                nTotalIonisationQReactions
                * cloud_.nParticle()
)/(counterIndex*deltaT*numberDensities_[0]
            * numberDensities_[1]*volume);

            Info<< "Ionisation reaction "
                <<  reactantMolA << " + " << reactantMolB
                <<  " --> "
                << reactantMolA << " + " << productMolE << " + " << productMolF
                << ", reaction rate = " << reactionRateIonisationQ
                << endl;
        }
    }
    else
    {
        label nTotalAssociativeIonisationReactions =
                    nTotalAssociativeIonisationReactions_;
        label nTotalIonisationPReactions = nTotalIonisationPReactions_;
        label nTotalIonisationQReactions = nTotalIonisationQReactions_;

        label nAssociativeIonisationReactionsPerTimeStep =
                    nAssociativeIonisationReactionsPerTimeStep_;
        label nIonisationPReactionsPerTimeStep =
                    nIonisationPReactionsPerTimeStep_;
        label nIonisationQReactionsPerTimeStep =
                    nIonisationQReactionsPerTimeStep_;

        if (Pstream::parRun())
        {
            reduce(nTotalAssociativeIonisationReactions, sumOp<label>());
            reduce(nTotalIonisationPReactions, sumOp<label>());
            reduce(nTotalIonisationQReactions, sumOp<label>());

            reduce(nAssociativeIonisationReactionsPerTimeStep, sumOp<label>());
            reduce(nIonisationPReactionsPerTimeStep, sumOp<label>());
            reduce(nIonisationQReactionsPerTimeStep, sumOp<label>());
        }

        if (nTotalAssociativeIonisationReactions > VSMALL)
        {
                word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
                word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];

                word productMolA =
                    cloud_.typeIdList()[associativeIonisationProductIds_[0]];
                word productMolB =
                    cloud_.typeIdList()[associativeIonisationProductIds_[1]];


                Info<< "Associative ionisation reaction "
                    <<  reactantMolA << " + " << reactantMolB
                    <<  " --> "
                    << productMolA << " + " << productMolB
                    << " is active, nReactions this time step = "
                    << nAssociativeIonisationReactionsPerTimeStep
                    << endl;
        }

        if (nTotalIonisationPReactions > VSMALL)
        {
                word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
                word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];

                word productMolA =
                        cloud_.typeIdList()[ionisationPProductIds_[0]];
                word productMolB =
                        cloud_.typeIdList()[ionisationPProductIds_[1]];


                Info<< "Ionisation reaction "
                    <<  reactantMolA << " + " << reactantMolB
                    <<  " --> "
                    << productMolA << " + " << productMolB << " + "
                    << reactantMolB
                    << " is active, nReactions this time step = "
                    << nIonisationPReactionsPerTimeStep
                    << endl;
        }

        if (nTotalIonisationQReactions > VSMALL)
        {
                word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
                word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];

                word productMolA =
                    cloud_.typeIdList()[ionisationQProductIds_[0]];
                word productMolB =
                    cloud_.typeIdList()[ionisationQProductIds_[1]];


                Info<< "Ionisation reaction "
                    <<  reactantMolA << " + " << reactantMolB
                    <<  " --> "
                    << reactantMolA << " + " << productMolA
                    << " + " << productMolB
                    << " is active, nReactions this time step = "
                    << nIonisationQReactionsPerTimeStep
                    << endl;
        }
    }

    nAssociativeIonisationReactionsPerTimeStep_ = 0.0;
    nIonisationPReactionsPerTimeStep_ = 0.0;
    nIonisationQReactionsPerTimeStep_ = 0.0;
}


bool forwardAssociativeIonisationDissimilarSpecies::relax() const
{
    return relax_;
}

} // End namespace Foam

// ************************************************************************* //
