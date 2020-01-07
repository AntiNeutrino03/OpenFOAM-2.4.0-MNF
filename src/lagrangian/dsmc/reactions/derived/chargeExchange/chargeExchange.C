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

#include "chargeExchange.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(chargeExchange, 0);

addToRunTimeSelectionTable(dsmcReaction, chargeExchange, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
chargeExchange::chargeExchange
(
    Time& t,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcReaction(t, cloud, dict),
    chargeExchangeProductIds_(),
    dissociationProductIds_(),
    ionisationProductIds_(),
    dissociationPossible_(Switch(propsDict_.lookup("dissociationPossible"))),
    heatOfReactionDissoc_(0.0),
    heatOfReactionIon_(
        readScalar(propsDict_.lookup("heatOfReactionIonisation"))),
    heatOfReactionExch_(
        readScalar(propsDict_.lookup("heatOfReactionExchange"))),
    activationEnergy_(0.0),
    aCoeff_(readScalar(propsDict_.lookup("aCoeff"))),
    bCoeff_(readScalar(propsDict_.lookup("bCoeff"))),
    nChargeExchangeReactions_(0),
    nIonisationReactions_(0),
    nDissociationReactions_(0),
    nChargeExchangeReactionsPerTimeStep_(0),
    nDissociationReactionsPerTimeStep_(0),
    nIonisationReactionsPerTimeStep_(0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

chargeExchange::~chargeExchange()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void chargeExchange::initialConfiguration()
{
    setCommonReactionProperties();
    setProperties();
}


void chargeExchange::setProperties()
{
    //check that reactants are 'ATOMS' or 'MOLECULES', not 'ELECTRONS'

    if (charge1_ == -1)
    {
        FatalErrorIn("chargeExchange::setProperties()")
            << "Reactant must not be an electron: " << reactants_[0]
            << nl
            << exit(FatalError);
    }

    if (vDof1_ > 1)
    {
            FatalErrorIn("chargeExchange::setProperties()")
        << "Reactions are currently only implemented for "
        << "monatomic and diatomic species"
        << " This is a polyatomic:" << reactants_[0]
        << nl
        << exit(FatalError);
    }
    
    if (charge2_ == -1)
    {
        FatalErrorIn("chargeExchange::setProperties()")
            << "Reactant must not be an electron: " << reactants_[1]
            << nl
            << exit(FatalError);
    }

    if (vDof2_ > 1)
    {
            FatalErrorIn("chargeExchange::setProperties()")
        << "Reactions are currently only implemented for "
        << "monatomic and diatomic species"
        << " This is a polyatomic:" << reactants_[1]
        << nl
        << exit(FatalError);
    }

    // reading in charge exchange products

    List<word> chargeExchangeProductMolecules
            (propsDict_.lookup("chargeExchangeProducts"));

    chargeExchangeProductIds_.setSize(
                            chargeExchangeProductMolecules.size(), -1);

    forAll(chargeExchangeProductIds_, i)
    {
        chargeExchangeProductIds_[i] = findIndex(cloud_.typeIdList(),
                                        chargeExchangeProductMolecules[i]);

        // check that products belong to the typeIdList
        if (chargeExchangeProductIds_[i] == -1)
        {
            FatalErrorIn("chargeExchange::setProperties()")
                << "Cannot find type id: "
                << chargeExchangeProductMolecules[i] << nl
                << exit(FatalError);
        }

        // check that products are a 'MOLECULE' or an 'ATOM'

        label charge =
                cloud_.constProps(chargeExchangeProductIds_[i]).charge();

        if (charge == -1)
        {
            FatalErrorIn("chargeExchange::setProperties()")
                << "Products cannot be an electron: "
                << chargeExchangeProductMolecules
                << nl
                << exit(FatalError);
        }
    }

    // reading in dissociation products

    if (dissociationPossible_)
    {
        const List<word> dissociationProductMolecules
                    (propsDict_.lookup("dissociationProducts"));

        dissociationProductIds_.setSize(
                                dissociationProductMolecules.size(), -1);

        forAll(dissociationProductIds_, i)
        {
            dissociationProductIds_[i] = findIndex(cloud_.typeIdList(),
                                            dissociationProductMolecules[i]);

            // check that products belong to the typeIdList
            if (dissociationProductIds_[i] == -1)
            {
                FatalErrorIn("chargeExchange::setProperties()")
                    << "Cannot find type id: "
                    << dissociationProductMolecules[i] << nl
                    << exit(FatalError);
            }

            // check that products are ATOMS

            label iD = dissociationProductIds_[i];

            scalar rDof =
                cloud_.constProps(iD).rotationalDegreesOfFreedom();

            if (rDof > VSMALL)
            {
                FatalErrorIn("chargeExchange::setProperties()")
                    << "Dissociation products must be atoms: "
                    << dissociationProductMolecules
                    << nl
                    << exit(FatalError);
            }
        }
    }

    // reading in ionisation products

    List<word> ionisationProductMolecules
                        (propsDict_.lookup("ionisationProducts"));

    ionisationProductIds_.setSize(ionisationProductMolecules.size(), -1);

    forAll(ionisationProductIds_, i)
    {
        ionisationProductIds_[i] = findIndex(cloud_.typeIdList(),
                                            ionisationProductMolecules[i]);

        // check that products belong to the typeIdList
        if (ionisationProductIds_[i] == -1)
        {
            FatalErrorIn("chargeExchange::setProperties()")
                << "Cannot find type id: "
                << ionisationProductMolecules[i] << nl
                << exit(FatalError);
        }
    }

    // check that second product is an electron

    const label charge = cloud_.constProps(ionisationProductIds_[1]).charge();

    if (charge != -1)
    {
        FatalErrorIn("chargeExchange::setProperties()")
            << "Second ionisation product must be an electron: "
            << ionisationProductMolecules
            << nl
            << exit(FatalError);
    }

    activationEnergy_ *= physicoChemical::k.value();

    if (dissociationPossible_)
    {
        heatOfReactionDissoc_ =
                readScalar(propsDict_.lookup("heatOfReactionDissociation"));
    }
}


bool chargeExchange::tryReactMolecules(const label typeIdP,
                                       const label typeIdQ)
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


void chargeExchange::reaction
(
    dsmcParcel& p,
    dsmcParcel& q,
    const DynamicList<label>& candidateList,
    const List<DynamicList<label> >& candidateSubList,
    const label candidateP,
    const List<label>& whichSubCell
)
{}


void chargeExchange::reaction
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

        bool dissocReactionP = false;
        bool ionisationReactionP = false;
        bool chargeExchange = false;

        //3 reactions possible
        // 1. Dissociation of P
        // 2. Ionisation of P
        // 3. Charge exchange


        scalar EcP = 0.0;
        scalar TColl = 0.0;

        // check for dissociation

        if (dissociationPossible_)
        {
            if (cloud_.constProps(typeIdP).rotationalDegreesOfFreedom() > 0)
            {
                label idP = reactionProperties(cloud_, p, q).charDissLevelP();
                label imaxP = 0;

                EcP = reactionProperties(cloud_, p, q).translationalEnergy() + reactionProperties(cloud_, p, q).EVibP();

                imaxP = EcP/(physicoChemical::k.value()*reactionProperties(cloud_, p, q).thetaVP());

                if (imaxP - idP > 0)
                {
                    // Dissociation can occur
                    totalReactionProbability += 1.0;
                    reactionProbabilities[0] = 1.0;
                }
            }
        }

        scalar ionisationEnergy =
                    cloud_.constProps(typeIdP).ionisationTemperature()
                    *physicoChemical::k.value();

        // calculate if an ionisation of species P is possible
        EcP = reactionProperties(cloud_, p, q).translationalEnergy() + reactionProperties(cloud_, p, q).EEleP();

        if ((EcP - ionisationEnergy) > VSMALL)
        {
            // Ionisation can occur
            totalReactionProbability += 1.0;
            reactionProbabilities[1] = 1.0;
        }

        scalar heatOfReactionExchJoules = heatOfReactionExch_
                                        *physicoChemical::k.value();

        scalar aDash = aCoeff_*(pow(2.5 - reactionProperties(cloud_, p, q).omegaPQ(), bCoeff_)
                        *exp(lgamma(2.5 - reactionProperties(cloud_, p, q).omegaPQ()))
                        /exp(lgamma(2.5 - reactionProperties(cloud_, p, q).omegaPQ() + bCoeff_)));

        TColl = (reactionProperties(cloud_, p, q).translationalEnergy()/(physicoChemical::k.value()))
                                        /(2.5 - reactionProperties(cloud_, p, q).omegaPQ());

        scalar activationEnergy = activationEnergy_
                            + (aDash*pow((TColl/273.0), bCoeff_)
                            * mag(heatOfReactionExchJoules));

         // scalar activationEnergy = activationEnergy_
                            //+ (aCoeff_*pow((5000/273.0), bCoeff_)
                           // *mag(heatOfReactionExchJoules));

        if (heatOfReactionExchJoules < 0.0)
        {
            activationEnergy -= heatOfReactionExchJoules;
        }

        if (EcP > activationEnergy)
        {
            label keyElectronicLevel = -1;

            for (label i = 0; i < (reactionProperties(cloud_, p, q).jMaxP() + 1); i++)
            {
                scalar electronicEnergy = reactionProperties(cloud_, p, q).EEListP()[i];

                if (electronicEnergy > activationEnergy)
                {
                    break;
                }

                keyElectronicLevel++;
            }

            label trialELevel = cloud_.postCollisionElectronicEnergyLevel
                            (
                                EcP,
                                (reactionProperties(cloud_, p, q).jMaxP() + 1),
                                reactionProperties(cloud_, p, q).omegaPQ(),
                                reactionProperties(cloud_, p, q).EEListP(),
                                reactionProperties(cloud_, p, q).gListP()
                            );

            if (trialELevel == keyElectronicLevel)
            {
                scalar prob = 0.0;

                label nPossStates = 0;

                if ((reactionProperties(cloud_, p, q).jMaxP() + 1) == 1)
                {
                    nPossStates = reactionProperties(cloud_, p, q).gListP()[0];
                }
                else
                {
                    forAll(reactionProperties(cloud_, p, q).EEListP(), n)
                    {
                        if (EcP > reactionProperties(cloud_, p, q).EEListP()[n])
                        {
                            nPossStates += reactionProperties(cloud_, p, q).gListP()[n];
                        }
                    }
                }

                label nState = ceil(cloud_.rndGen().scalar01()*(nPossStates));
                label nAvailableStates = 0;
                label nLevel = -1;

                forAll(reactionProperties(cloud_, p, q).EEListP(), n)
                {
                    nAvailableStates += reactionProperties(cloud_, p, q).gListP()[n];

                    if (nState <= nAvailableStates && nLevel < 0)
                    {
                        nLevel = n;
                    }
                }

                // Calculate the probability of it occurring
                scalar summation = 0.0;

                for (label i = 0; i <= nLevel; i++)
                {
                    summation += reactionProperties(cloud_, p, q).gListP()[i]*pow((EcP - reactionProperties(cloud_, p, q).EEListP()[i]), 1.5 - reactionProperties(cloud_, p, q).omegaPQ());
                }

                prob = (reactionProperties(cloud_, p, q).gListP()[trialELevel]
                    *pow((EcP - reactionProperties(cloud_, p, q).EEListP()[trialELevel]), 1.5 - reactionProperties(cloud_, p, q).omegaPQ()))/summation;

                if (prob > cloud_.rndGen().scalar01())
                {
                    // Charge exchange can occur
                    totalReactionProbability += prob;
                    reactionProbabilities[2] = prob;
                }
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
                        if (i == 2)
                        {
                            // Dissociation is to occur
                            chargeExchange = true;
                            break;
                        }
                    }
                }
            }
        }

        if (dissocReactionP)
        {
            nDissociationReactions_++;
            nDissociationReactionsPerTimeStep_++;

            if (allowSplitting_)
            {
                relax_ = false;
                
                reactionProperties(cloud_, p, q).dissociateP
                (
                    heatOfReactionDissoc_,
                    dissociationProductIds_,
                    p,
                    q
                );
            }
        }

        if (ionisationReactionP)
        {
            nIonisationReactions_++;
            nIonisationReactionsPerTimeStep_++;

            if (allowSplitting_)
            {
                relax_ = false;
                
                reactionProperties(cloud_, p, q).ioniseP
                (
                    heatOfReactionIon_,
                    ionisationProductIds_,
                    p,
                    q
                );
            }
        }

        if (chargeExchange)
        {
            nChargeExchangeReactions_++;
            nChargeExchangeReactionsPerTimeStep_++;

            if (allowSplitting_)
            {
                relax_ = false;
                
                scalar translationalEnergy = reactionProperties(cloud_, p, q).translationalEnergy();

                translationalEnergy += heatOfReactionExchJoules;

                translationalEnergy += reactionProperties(cloud_, p, q).ERotP() + reactionProperties(cloud_, p, q).EVibP() + reactionProperties(cloud_, p, q).EEleP()
                                        + reactionProperties(cloud_, p, q).ERotQ() + reactionProperties(cloud_, p, q).EVibQ() + reactionProperties(cloud_, p, q).EEleQ();

                scalar mR = cloud_.constProps(chargeExchangeProductIds_[0]).mass()*cloud_.constProps(chargeExchangeProductIds_[1]).mass()/(cloud_.constProps(chargeExchangeProductIds_[0]).mass() + cloud_.constProps(chargeExchangeProductIds_[1]).mass());

                scalar relVel = sqrt((2.0*translationalEnergy)/mR);

                // Variable Hard Sphere collision part for collision of
                // molecules
                scalar cosTheta = 2.0*cloud_.rndGen().scalar01() - 1.0;

                scalar sinTheta = sqrt(1.0 - cosTheta*cosTheta);

                scalar phi = twoPi*cloud_.rndGen().scalar01();

                vector postCollisionRelU =
                    relVel
                    *vector
                        (
                            cosTheta,
                            sinTheta*cos(phi),
                            sinTheta*sin(phi)
                        );

                scalar mP = cloud_.constProps(chargeExchangeProductIds_[0]).mass();
                scalar mQ = cloud_.constProps(chargeExchangeProductIds_[1]).mass();

                vector UP = reactionProperties(cloud_, p, q).Ucm() + (postCollisionRelU*mQ/(mP + mQ));
                vector UQ = reactionProperties(cloud_, p, q).Ucm() - (postCollisionRelU*mP/(mP + mQ));
                
                vector positionP(p.position());
    
                label cell = -1;
                label tetFace = -1;
                label tetPt = -1;

                mesh_.findCellFacePt
                (
                    positionP,
                    cell,
                    tetFace,
                    tetPt
                );
                
                scalar RWFp = p.RWF();
                scalar ERotP = 0.0;
                labelList vibLevelP(0, 0);

                if (cloud_.constProps(chargeExchangeProductIds_[0])
                    .vibrationalDegreesOfFreedom()
                > VSMALL)
                {
                    vibLevelP.setSize(1);
                    vibLevelP[0] = 0;
                }
                label ELevelP = 0;
                
                p.position() = positionP;
                p.U() = UP;
                p.RWF() = RWFp;
                p.ERot() = ERotP;
                p.ELevel() = ELevelP;
                p.cell() = cell;
                p.tetFace() = tetFace;
                p.tetPt() = tetPt;
                p.typeId() = chargeExchangeProductIds_[0];
                p.newParcel() = 0;
                p.vibLevel() = vibLevelP;
                
                vector positionQ(q.position());
    
                cell = -1;
                tetFace = -1;
                tetPt = -1;

                mesh_.findCellFacePt
                (
                    positionQ,
                    cell,
                    tetFace,
                    tetPt
                );
                
                scalar RWFq = q.RWF();
                scalar ERotQ = 0.0;
                labelList vibLevelQ(0, 0);

                if (cloud_.constProps(chargeExchangeProductIds_[1])
                    .vibrationalDegreesOfFreedom()
                > VSMALL)
                {
                    vibLevelQ.setSize(1);
                    vibLevelQ[0] = 0;
                }
                label ELevelQ = 0;
                
                q.position() = positionQ;
                q.U() = UQ;
                q.RWF() = RWFq;
                q.ERot() = ERotQ;
                q.ELevel() = ELevelQ;
                q.cell() = cell;
                q.tetFace() = tetFace;
                q.tetPt() = tetPt;
                q.typeId() = chargeExchangeProductIds_[1];
                q.newParcel() = 0;
                q.vibLevel() = vibLevelQ;
            }
        }
    }

    if (typeIdP == reactantIds_[1] && typeIdQ == reactantIds_[0])
    {
        relax_ = true;

        scalar totalReactionProbability = 0.0;
        scalarList reactionProbabilities(3, 0.0);

        bool dissocReactionQ = false;
        bool ionisationReactionQ = false;
        bool chargeExchange = false;

        //3 reactions possible
        // 1. Dissociation of P
        // 2. Ionisation of P
        // 3. Charge exchange


        scalar EcQ = 0.0;
        scalar TColl = 0.0;

        // check for dissociation

        if (dissociationPossible_)
        {
            if (cloud_.constProps(typeIdQ).rotationalDegreesOfFreedom() > 0)
            {
                label idQ = reactionProperties(cloud_, p, q).charDissLevelQ();
                label imaxQ = 0;

                EcQ = reactionProperties(cloud_, p, q).translationalEnergy() + reactionProperties(cloud_, p, q).EVibQ();

                imaxQ = EcQ/(physicoChemical::k.value()*reactionProperties(cloud_, p, q).thetaVQ());

                if (imaxQ - idQ > 0)
                {
                    // Dissociation can occur
                    totalReactionProbability += 1.0;
                    reactionProbabilities[0] = 1.0;
                }
            }
        }

        scalar ionisationEnergy =
                    cloud_.constProps(typeIdQ).ionisationTemperature()
                    *physicoChemical::k.value();

        // calculate if an ionisation of species Q is possible
        EcQ = reactionProperties(cloud_, p, q).translationalEnergy() + reactionProperties(cloud_, p, q).EEleQ();

        if ((EcQ - ionisationEnergy) > VSMALL)
        {
            // Ionisation can occur
            totalReactionProbability += 1.0;
            reactionProbabilities[1] = 1.0;
        }

        scalar heatOfReactionExchJoules = heatOfReactionExch_
                *physicoChemical::k.value();

        scalar aDash = aCoeff_*(pow(2.5 - reactionProperties(cloud_, p, q).omegaPQ(), bCoeff_)
                        *exp(lgamma(2.5 - reactionProperties(cloud_, p, q).omegaPQ()))
                        /exp(lgamma(2.5 - reactionProperties(cloud_, p, q).omegaPQ() + bCoeff_)));

        TColl = (reactionProperties(cloud_, p, q).translationalEnergy()/(physicoChemical::k.value()))
                                            /(2.5 - reactionProperties(cloud_, p, q).omegaPQ());

        scalar activationEnergy = (aDash*pow((TColl/273.0), bCoeff_)
                    * mag(heatOfReactionExchJoules));

//         scalar activationEnergy = (aCoeff_*pow((5000/273.0), bCoeff_)
                            //* mag(heatOfReactionExchJoules));

        if (heatOfReactionExchJoules < 0.0)
        {
            activationEnergy -= heatOfReactionExchJoules;
        }

        if (EcQ > activationEnergy)
        {
            label keyElectronicLevel = -1;

            for (label i = 0; i < (reactionProperties(cloud_, p, q).jMaxQ() + 1); i++)
            {
                scalar electronicEnergy = reactionProperties(cloud_, p, q).EEListQ()[i];

                if (electronicEnergy > activationEnergy)
                {
                    break;
                }

                keyElectronicLevel++;
            }

            label trialELevel = cloud_.postCollisionElectronicEnergyLevel
                            (
                                EcQ,
                                (reactionProperties(cloud_, p, q).jMaxQ() + 1),
                                reactionProperties(cloud_, p, q).omegaPQ(),
                                reactionProperties(cloud_, p, q).EEListQ(),
                                reactionProperties(cloud_, p, q).gListQ()
                            );

            if (trialELevel == keyElectronicLevel)
            {
                scalar prob = 0.0;

                label nPossStates = 0;

                if ((reactionProperties(cloud_, p, q).jMaxQ() + 1) == 1)
                {
                    nPossStates = reactionProperties(cloud_, p, q).gListQ()[0];
                }
                else
                {
                    forAll(reactionProperties(cloud_, p, q).EEListQ(), n)
                    {
                        if (EcQ > reactionProperties(cloud_, p, q).EEListQ()[n])
                        {
                            nPossStates += reactionProperties(cloud_, p, q).gListQ()[n];
                        }
                    }
                }

                label nState = ceil(cloud_.rndGen().scalar01()*(nPossStates));
                label nAvailableStates = 0;
                label nLevel = -1;

                forAll(reactionProperties(cloud_, p, q).EEListQ(), n)
                {
                    nAvailableStates += reactionProperties(cloud_, p, q).gListQ()[n];

                    if (nState <= nAvailableStates && nLevel < 0)
                    {
                        nLevel = n;
                    }
                }

                // Calculate the probability of it occurring

                scalar summation = 0.0;

                for (label i = 0; i <= nLevel; i++)
                {
                    summation += reactionProperties(cloud_, p, q).gListQ()[i]*pow((EcQ - reactionProperties(cloud_, p, q).EEListQ()[i]), 1.5 - reactionProperties(cloud_, p, q).omegaPQ());
                }

                prob = (reactionProperties(cloud_, p, q).gListQ()[trialELevel]
                    *pow((EcQ - reactionProperties(cloud_, p, q).EEListQ()[trialELevel]), 1.5 - reactionProperties(cloud_, p, q).omegaPQ()))/summation;

                if (prob > cloud_.rndGen().scalar01())
                {
                    // Charge exchange can occur
                    totalReactionProbability += prob;
                    reactionProbabilities[2] = prob;
                }
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
                            // Dissociation is to occur
                            dissocReactionQ = true;
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
                            // Dissociation is to occur
                            chargeExchange = true;
                            break;
                        }
                    }
                }
            }
        }

        if (dissocReactionQ)
        {
            nDissociationReactions_++;
            nDissociationReactionsPerTimeStep_++;

            if (allowSplitting_)
            {
                relax_ = false;

                reactionProperties(cloud_, p, q).dissociateQ
                (
                    heatOfReactionDissoc_,
                    dissociationProductIds_,
                    p,
                    q
                );
            }
        }

        if (ionisationReactionQ)
        {
            nIonisationReactions_++;
            nIonisationReactionsPerTimeStep_++;

            if (allowSplitting_)
            {
                relax_ = false;

                reactionProperties(cloud_, p, q).ioniseQ
                (
                    heatOfReactionIon_,
                    ionisationProductIds_,
                    p,
                    q
                );
            }
        }

        if (chargeExchange)
        {
            nChargeExchangeReactions_++;
            nChargeExchangeReactionsPerTimeStep_++;

            if (allowSplitting_)
            {
                relax_ = false;

                scalar translationalEnergy = reactionProperties(cloud_, p, q).translationalEnergy();

                translationalEnergy += heatOfReactionExchJoules;

                translationalEnergy += reactionProperties(cloud_, p, q).ERotP() + reactionProperties(cloud_, p, q).EVibP() + reactionProperties(cloud_, p, q).EEleP()
                                    + reactionProperties(cloud_, p, q).ERotQ() + reactionProperties(cloud_, p, q).EVibQ() + reactionProperties(cloud_, p, q).EEleQ();
                                    
                scalar mR = cloud_.constProps(chargeExchangeProductIds_[0]).mass()*cloud_.constProps(chargeExchangeProductIds_[1]).mass()/(cloud_.constProps(chargeExchangeProductIds_[0]).mass() + cloud_.constProps(chargeExchangeProductIds_[1]).mass());

                scalar relVel = sqrt((2.0*translationalEnergy)/mR);

                // Variable Hard Sphere collision part for collision of
                // molecules
                scalar cosTheta = 2.0*cloud_.rndGen().scalar01() - 1.0;

                scalar sinTheta = sqrt(1.0 - cosTheta*cosTheta);

                scalar phi = twoPi*cloud_.rndGen().scalar01();

                vector postCollisionRelU =
                    relVel
                    *vector
                        (
                            cosTheta,
                            sinTheta*cos(phi),
                            sinTheta*sin(phi)
                        );

                scalar mP = cloud_.constProps(chargeExchangeProductIds_[1]).mass();
                scalar mQ = cloud_.constProps(chargeExchangeProductIds_[0]).mass();

                vector UP = reactionProperties(cloud_, p, q).Ucm() + (postCollisionRelU*mQ/(mP + mQ));
                vector UQ = reactionProperties(cloud_, p, q).Ucm() - (postCollisionRelU*mP/(mP + mQ));
                
                vector positionP(p.position());
    
                label cell = -1;
                label tetFace = -1;
                label tetPt = -1;

                mesh_.findCellFacePt
                (
                    positionP,
                    cell,
                    tetFace,
                    tetPt
                );
                
                scalar RWFp = p.RWF();
                scalar ERotP = 0.0;
                labelList vibLevelP(0, 0);

                if (cloud_.constProps(chargeExchangeProductIds_[1])
                    .vibrationalDegreesOfFreedom()
                > VSMALL)
                {
                    vibLevelP.setSize(1);
                    vibLevelP[0] = 0;
                }
                label ELevelP = 0;
                
                p.position() = positionP;
                p.U() = UP;
                p.RWF() = RWFp;
                p.ERot() = ERotP;
                p.ELevel() = ELevelP;
                p.cell() = cell;
                p.tetFace() = tetFace;
                p.tetPt() = tetPt;
                p.typeId() = chargeExchangeProductIds_[1];
                p.newParcel() = 0;
                p.vibLevel() = vibLevelP;
                
                vector positionQ(q.position());
    
                cell = -1;
                tetFace = -1;
                tetPt = -1;

                mesh_.findCellFacePt
                (
                    positionQ,
                    cell,
                    tetFace,
                    tetPt
                );
                
                scalar RWFq = q.RWF();
                scalar ERotQ = 0.0;
                labelList vibLevelQ(0, 0);

                if (cloud_.constProps(chargeExchangeProductIds_[0])
                    .vibrationalDegreesOfFreedom()
                > VSMALL)
                {
                    vibLevelQ.setSize(1);
                    vibLevelQ[0] = 0;
                }
                label ELevelQ = 0;
                
                q.position() = positionQ;
                q.U() = UQ;
                q.RWF() = RWFq;
                q.ERot() = ERotQ;
                q.ELevel() = ELevelQ;
                q.cell() = cell;
                q.tetFace() = tetFace;
                q.tetPt() = tetPt;
                q.typeId() = chargeExchangeProductIds_[0];
                q.newParcel() = 0;
                q.vibLevel() = vibLevelQ;
            }
        }
    }
}


void  chargeExchange::outputResults(const label counterIndex)
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
        label nTotChargeExchangeReactions = nChargeExchangeReactions_;
        label nTotDissociationReactions = nDissociationReactions_;
        label nTotIonisationReactions = nIonisationReactions_;

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
            reduce(nTotChargeExchangeReactions, sumOp<label>());
            reduce(nTotDissociationReactions, sumOp<label>());
            reduce(nTotIonisationReactions, sumOp<label>());
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

        word productMolA = cloud_.typeIdList()[chargeExchangeProductIds_[0]];
        word productMolB = cloud_.typeIdList()[chargeExchangeProductIds_[1]];

        word productMolC;
        word productMolD;

        if (dissociationPossible_)
        {
            productMolC = cloud_.typeIdList()[dissociationProductIds_[0]];
            productMolD = cloud_.typeIdList()[dissociationProductIds_[1]];
        }

        word productMolE = cloud_.typeIdList()[ionisationProductIds_[0]];
        word productMolF = cloud_.typeIdList()[ionisationProductIds_[1]];

        if ((numberDensities_[0] > 0.0) && (numberDensities_[1] > 0.0))
        {
            scalar reactionRateChargeExchange = 0.0;
            scalar reactionRateDissociation = 0.0;
            scalar reactionRateIonisation = 0.0;

            reactionRateChargeExchange =
            (
                nTotChargeExchangeReactions
                * cloud_.nParticle()
)/(counterIndex*deltaT*numberDensities_[0]
            * numberDensities_[1]*volume);

            Info<< "Charge exchange reaction "
                <<  reactantMolA << " + " << reactantMolB
                <<  " --> "
                << productMolA << " + " << productMolB
                << ", reaction rate = " << reactionRateChargeExchange
                << endl;

            reactionRateIonisation =
            (
                nTotIonisationReactions
                * cloud_.nParticle()
)/(counterIndex*deltaT*numberDensities_[0]
            * numberDensities_[1]*volume);

            Info<< "Ionisation reaction "
                <<  reactantMolA << " + " << reactantMolB
                <<  " --> "
                << productMolE << " + " << productMolF
                 << " + " << reactantMolB 
                << ", reaction rate = " << reactionRateIonisation
                << endl;

            if (dissociationPossible_)
            {
                reactionRateDissociation =
                (
                    nTotDissociationReactions
                    * cloud_.nParticle()
)/(counterIndex*deltaT*numberDensities_[0]
                * numberDensities_[1]*volume);

                Info<< "Dissociation reaction "
                    <<  reactantMolA << " + " << reactantMolB
                    <<  " --> "
                    << productMolC << " + " << productMolD
                    << " + " << reactantMolB 
                    << ", reaction rate = " << reactionRateDissociation
                    << endl;
            }
        }
    }
    else
    {
        label nTotChargeExchangeReactions = nChargeExchangeReactions_;
        label nTotDissociationReactions = nDissociationReactions_;
        label nTotIonisationReactions = nIonisationReactions_;

        label nChargeExchangeReactionsPerTimeStep =
                    nChargeExchangeReactionsPerTimeStep_;
        label nDissociationReactionsPerTimeStep =
                    nDissociationReactionsPerTimeStep_;
        label nIonisationReactionsPerTimeStep =
                    nIonisationReactionsPerTimeStep_;


        if (Pstream::parRun())
        {
            reduce(nTotChargeExchangeReactions, sumOp<label>());
            reduce(nTotDissociationReactions, sumOp<label>());
            reduce(nTotIonisationReactions, sumOp<label>());
            reduce(nChargeExchangeReactionsPerTimeStep, sumOp<label>());
            reduce(nDissociationReactionsPerTimeStep, sumOp<label>());
            reduce(nIonisationReactionsPerTimeStep, sumOp<label>());
        }

        if (nTotChargeExchangeReactions > VSMALL)
        {
                word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
                word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];

                word productMolA =
                        cloud_.typeIdList()[chargeExchangeProductIds_[0]];
                word productMolB =
                        cloud_.typeIdList()[chargeExchangeProductIds_[1]];

                Info<< "Charge exchange reaction "
                    <<  reactantMolA << " + " << reactantMolB
                    <<  " --> "
                    << productMolA << " + " << productMolB
                    << " is active, nReactions this time step = "
                    << nChargeExchangeReactionsPerTimeStep << endl;
        }

        if (nTotDissociationReactions > VSMALL)
        {
                word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
                word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];

                word productMolA =
                            cloud_.typeIdList()[dissociationProductIds_[0]];
                word productMolB =
                            cloud_.typeIdList()[dissociationProductIds_[1]];

                Info<< "Dissociation reaction "
                    <<  reactantMolA << " + " << reactantMolB
                    <<  " --> "
                    << productMolA << " + " << productMolB
                    << " + " << reactantMolB
                    << " is active, nReactions this time step = "
                    << nDissociationReactionsPerTimeStep << endl;
        }

        if (nTotIonisationReactions > VSMALL)
        {
                word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
                word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];

                word productMolA =
                        cloud_.typeIdList()[ionisationProductIds_[0]];
                word productMolB =
                        cloud_.typeIdList()[ionisationProductIds_[1]];

                Info<< "Ionisation reaction "
                    <<  reactantMolA << " + " << reactantMolB
                    <<  " --> "
                    << productMolA << " + " << productMolB
                    << " + " << reactantMolB
                    << " is active, nReactions this time step = "
                    << nIonisationReactionsPerTimeStep << endl;
        }
    }

    nChargeExchangeReactionsPerTimeStep_ = 0.0;
    nDissociationReactionsPerTimeStep_ = 0.0;
    nIonisationReactionsPerTimeStep_ = 0.0;
}


bool chargeExchange::relax() const
{
    return relax_;
}

} // End namespace Foam

// ************************************************************************* //
