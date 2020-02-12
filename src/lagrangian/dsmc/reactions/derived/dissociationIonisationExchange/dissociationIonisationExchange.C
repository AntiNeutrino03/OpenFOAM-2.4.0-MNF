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

#include "dissociationIonisationExchange.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(dissociationIonisationExchange, 0);

addToRunTimeSelectionTable
(
    dsmcReaction,
    dissociationIonisationExchange,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

dissociationIonisationExchange::dissociationIonisationExchange
(
    Time& t,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcReaction(t, cloud, dict),
    exchangeProductIds_(),
    chargeExchangeProductIds_(),
    dissociationProductIds_(),
    ionisationProductsIdsP_(),
    ionisationProductsIdsQ_(),
    chargedAtom_(Switch(propsDict_.lookup("chargedAtom"))),
    chargedMolecule_(Switch(propsDict_.lookup("chargedMolecule"))),
    chargeExchange_(Switch(propsDict_.lookup("chargeExchange"))),
    heatOfReactionDiss_(),
    heatOfReactionExch_(readScalar(propsDict_.lookup("heatOfReactionExch"))),
    heatOfReactionChargeExchange_(),
    heatOfReactionIonP_(),
    heatOfReactionIonQ_(),
    aCoeff_(readScalar(propsDict_.lookup("aCoeff"))),
    bCoeff_(readScalar(propsDict_.lookup("bCoeff"))),
    aCoeffCharge_(),
    bCoeffCharge_(),
    nTotExchangeReactions_(0),
    nTotChargeExchangeReactions_(),
    nTotDissociationReactions_(0),
    nTotIonisationReactionsP_(0),
    nTotIonisationReactionsQ_(0),
    nExchangeReactionsPerTimeStep_(0),
    nChargeExchangeReactionsPerTimeStep_(0),
    nDissociationReactionsPerTimeStep_(0),
    nIonisationReactionsPPerTimeStep_(0),
    nIonisationReactionsQPerTimeStep_(0)
{
    if (!chargedMolecule_)
    {
        heatOfReactionDiss_ = readScalar(
            propsDict_.lookup("heatOfReactionDiss"));
        heatOfReactionIonP_ = readScalar(
            propsDict_.lookup("heatOfReactionIonP"));
    }

    if (!chargedAtom_)
    {
        heatOfReactionIonQ_ = readScalar(
            propsDict_.lookup("heatOfReactionIonQ"));
    }

    if (chargeExchange_)
    {
        heatOfReactionChargeExchange_ = readScalar(
            propsDict_.lookup("heatOfReactionChargeExch"));
        aCoeffCharge_ = readScalar(
            propsDict_.lookup("aCoeffChargeExchange"));
        bCoeffCharge_ = readScalar(
            propsDict_.lookup("bCoeffChargeExchange"));
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dissociationIonisationExchange::~dissociationIonisationExchange()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dissociationIonisationExchange::initialConfiguration()
{
    setCommonReactionProperties();
    setProperties();
}


void dissociationIonisationExchange::setProperties()
{
    // check that the first reactant is a 'MOLECULE'

    if (rDof1_ < VSMALL)
    {
        FatalErrorIn("dissociationIonisationExchange::setProperties()")
            << "Reactant 1 must be a molecule (not an atom): "
            << reactants_[0]
            << nl
            << exit(FatalError);
    }

    if (vDof1_ > 1)
    {
         FatalErrorIn("dissociationIonisationExchange::setProperties()")
            << "Reactions are currently only implemented for "
            << "monatomic and diatomic species."
            << " This is a polyatomic:" << reactants_[0]
            << nl
            << exit(FatalError);
    }

    // check that the second reactant is an 'ATOM'

    if (rDof2_ > VSMALL)
    {
        FatalErrorIn("dissociationIonisationExchange::setProperties()")
            << "Reactant 2 must be an atom (not a molecule): "
            << reactants_[1]
            << nl
            << exit(FatalError);
    }

    //reading in products

    const List<word> exchangeProductMolecules
    (propsDict_.lookup("productsOfExchangeReaction"));

    if (exchangeProductMolecules.size() != 2)
    {
        FatalErrorIn("dissociationIonisationExchange::setProperties()")
            << "There should be two dissociationIonisationExchange "
            << "reaction products, instead of "
            << reactants_.size() << nl
            << exit(FatalError);
    }

    if (exchangeProductMolecules[0] == exchangeProductMolecules[1])
    {
        FatalErrorIn("dissociationIonisationExchange::setProperties()")
            << "Exchange reaction product molecules cannot be same species."
            << nl
            << exit(FatalError);
    }

    exchangeProductIds_.setSize(exchangeProductMolecules.size(), -1);

    forAll(exchangeProductIds_, r)
    {
        exchangeProductIds_[r] = findIndex(cloud_.typeIdList(),
                                           exchangeProductMolecules[r]);

        // check that reactants belong to the typeIdList
        if (exchangeProductIds_[r] == -1)
        {
            FatalErrorIn("dissociationIonisationExchange::setProperties()")
                << "Cannot find type id: " << exchangeProductMolecules[r]
                << nl
                << exit(FatalError);
        }
    }

    // check that first exchange product is a 'MOLECULE' (not an 'ATOM')

    const scalar rDof3 =
        cloud_.constProps(exchangeProductIds_[0]).rotationalDegreesOfFreedom();
        
    const scalar vDof3 =
        cloud_.constProps(exchangeProductIds_[0]).vibrationalDegreesOfFreedom();

    if (rDof3 < 0)
    {
        FatalErrorIn("dissociationIonisationExchange::setProperties()")
            << "First product of the exchange reaction must"
            << " be a molecule (not an atom): "
            << exchangeProductMolecules[0]
            << nl
            << exit(FatalError);
    }
    
    if (vDof3 > 1)
    {
         FatalErrorIn("dissociationIonisationExchange::setProperties()")
            << "Reactions are currently only implemented for "
            << "monatomic and diatomic species."
            << " This is a polyatomic:" << exchangeProductMolecules[0]
            << nl
            << exit(FatalError);
    }

    // check that second exchange product is an 'ATOM' (not a 'MOLECULE')

    label id = exchangeProductIds_[1];

    const scalar rDof4 = cloud_.constProps(id).rotationalDegreesOfFreedom();

    if (rDof4 > 0)
    {
        FatalErrorIn("dissociationIonisationExchange::setProperties()")
            << "Second product of the exchange reaction "
            << "must be an atom (not a molecule): "
            << exchangeProductMolecules[1]
            << nl
            << exit(FatalError);
    }

    if (chargeExchange_)
    {
        // reading in products

        const List<word> chargeExchangeProductMolecules
        (propsDict_.lookup("productsOfChargeExchange"));

        if (chargeExchangeProductMolecules.size() != 2)
        {
            FatalErrorIn("dissociationIonisationExchange::setProperties()")
                << "There should be two charge exchange "
                << "reaction products, instead of "
                << chargeExchangeProductMolecules.size() << nl
                << exit(FatalError);
        }

        if (chargeExchangeProductMolecules[0]
            == chargeExchangeProductMolecules[1])
        {
            FatalErrorIn("dissociationIonisationExchange::setProperties()")
                << "Charge exchange reaction product molecules "
                << "cannot be same species." << nl
                << exit(FatalError);
        }

        chargeExchangeProductIds_.setSize
        (chargeExchangeProductMolecules.size(), -1);

        forAll(chargeExchangeProductIds_, r)
        {
            chargeExchangeProductIds_[r] = findIndex(cloud_.typeIdList(),
                                            chargeExchangeProductMolecules[r]);

            // check that reactants belong to the typeIdList
            if (chargeExchangeProductIds_[r] == -1)
            {
                FatalErrorIn("dissociationIonisationExchange::setProperties()")
                    << "Cannot find type id: "
                    << chargeExchangeProductMolecules[r]
                    << nl
                    << exit(FatalError);
            }
        }

        // check that first exchange product is a 'MOLECULE'

        label id2 = chargeExchangeProductIds_[0];

        const scalar rDof5 =
                    cloud_.constProps(id2).rotationalDegreesOfFreedom();

        if (rDof5 < VSMALL)
        {
            FatalErrorIn("dissociationIonisationExchange::setProperties()")
                << "First product of the exchange reaction must "
                << "be a molecule (not an atom): "
                << chargeExchangeProductMolecules[0]
                << nl
                << exit(FatalError);
        }

        // check that second exchange product is an 'ATOM'

        label id3 = chargeExchangeProductIds_[1];

        const scalar rDof6 =
                    cloud_.constProps(id3).rotationalDegreesOfFreedom();

        if (rDof6 > VSMALL)
        {
            FatalErrorIn("dissociationIonisationExchange::setProperties()")
                << "Second product of the exchange reaction "
                << "must be an atom (not a molecule): "
                << chargeExchangeProductMolecules[1]
                << nl
                << exit(FatalError);
        }
    }

    if (!chargedMolecule_)
    {
        const List<word> dissociationProductMolecules
        (propsDict_.lookup("productsOfDissociatedMolecule"));

        if (dissociationProductMolecules.size() != 2)
        {
            FatalErrorIn("dissociationIonisationExchange::setProperties()")
                << "There should be two dissociation reaction products, "
                << "instead of "
                << dissociationProductMolecules.size() << nl
                << exit(FatalError);
        }

        dissociationProductIds_.setSize
        (dissociationProductMolecules.size(), -1);

        forAll(dissociationProductIds_, r)
        {
            dissociationProductIds_[r] = findIndex(cloud_.typeIdList(),
                                             dissociationProductMolecules[r]);

            // check that reactants belong to the typeIdList
            if (dissociationProductIds_[r] == -1)
            {
                FatalErrorIn("dissociationIonisationExchange::setProperties()")
                    << "Cannot find type id: "
                    << dissociationProductMolecules[r]
                    << nl
                    << exit(FatalError);
            }
        }

        label id4 = dissociationProductIds_[0];

        const scalar rDof7=
                cloud_.constProps(id4).rotationalDegreesOfFreedom();

        if (rDof7 > VSMALL)
        {
            FatalErrorIn("dissociationIonisationExchange::setProperties()")
                << "First product of the dissociation reaction must "
                << "be an atom (not a molecule): "
                << dissociationProductMolecules[0]
                << nl
                << exit(FatalError);
        }

        // check that second exchange product is an 'ATOM' (not a 'MOLECULE')

        label id5 = dissociationProductIds_[1];

        const scalar rDof8 =
                cloud_.constProps(id5).rotationalDegreesOfFreedom();

        if (rDof8 > VSMALL)
        {
            FatalErrorIn("dissociationIonisationExchange::setProperties()")
                << "Second product of the exchange reaction must "
                << "be an atom (not a molecule): "
                << dissociationProductMolecules[1]
                << nl
                << exit(FatalError);
        }

        // read in ionisation products

        const List<word> ionisationProductMoleculesP
        (propsDict_.lookup("productsOfIonisedMolecule"));

        if (ionisationProductMoleculesP.size() != 2)
        {
            FatalErrorIn("dissociationIonisationExchange::setProperties()")
                << "There should be two ionisation reaction products, "
                << "instead of "
                << ionisationProductMoleculesP.size() << nl
                << exit(FatalError);
        }

        ionisationProductsIdsP_.setSize
        (ionisationProductMoleculesP.size(), -1);

        forAll(ionisationProductsIdsP_, r)
        {
            ionisationProductsIdsP_[r] = findIndex(cloud_.typeIdList(),
                                               ionisationProductMoleculesP[r]);

            // check that reactants belong to the typeIdList
            if (ionisationProductsIdsP_[r] == -1)
            {
                FatalErrorIn("dissociationIonisationExchange::setProperties()")
                    << "Cannot find type id: "
                    << ionisationProductMoleculesP[r]
                    << nl
                    << exit(FatalError);
            }
        }

        // check that first ionisation product is a 'MOLECULE' (not an 'ATOM')

        label id6 = ionisationProductsIdsP_[0];

        const scalar rDof9 =
                cloud_.constProps(id6).rotationalDegreesOfFreedom();

        if (rDof9 < VSMALL)
        {
            FatalErrorIn("dissociationIonisationExchange::setProperties()")
                << "First product of the molecule ionisation reaction "
                << "must be a molecule (not an atom): "
                << ionisationProductMoleculesP[0]
                << nl
                << exit(FatalError);
        }

        // check that second ionisation product is a 'ELECTRON'

        label id7 = ionisationProductsIdsP_[1];

        const label charge = cloud_.constProps(id7).charge();

        if (charge != -1)
        {
            FatalErrorIn("dissociationIonisationExchange::setProperties()")
                << "Second product of the molecule ionisation reaction must "
                << "be an electron: " << ionisationProductMoleculesP[1]
                << nl
                << exit(FatalError);
        }
    }

    if (!chargedAtom_)
    {
        const List<word> ionisationProductMoleculesQ
        (propsDict_.lookup("productsOfIonisedAtom"));

        if (ionisationProductMoleculesQ.size() != 2)
        {
            FatalErrorIn("dissociationIonisationExchange::setProperties()")
                << "There should be two ionisation reaction products, "
                << "instead of "
                << ionisationProductMoleculesQ.size() << nl
                << exit(FatalError);
        }

        ionisationProductsIdsQ_.setSize
        (ionisationProductMoleculesQ.size(), -1);

        forAll(ionisationProductsIdsQ_, r)
        {
            ionisationProductsIdsQ_[r] = findIndex(cloud_.typeIdList(),
                                               ionisationProductMoleculesQ[r]);

            // check that reactants belong to the typeIdList
            if (ionisationProductsIdsQ_[r] == -1)
            {
                FatalErrorIn("dissociationIonisationExchange::setProperties()")
                    << "Cannot find type id: "
                    << ionisationProductMoleculesQ[r]
                    << nl
                    << exit(FatalError);
            }
        }

        // check that first ionisation product is an 'ATOM'

        label id8 = ionisationProductsIdsQ_[0];

        const scalar rDof10 =
                cloud_.constProps(id8).rotationalDegreesOfFreedom();

        if (rDof10 > VSMALL)
        {
            FatalErrorIn("dissociationIonisationExchange::setProperties()")
                << "First product of the atom ionisation reaction must be "
                << "an atom (not a molecule): "
                << ionisationProductMoleculesQ[0]
                << nl
                << exit(FatalError);
        }

        // check that second ionisation product is a 'ELECTRON'

        label id9 = ionisationProductsIdsQ_[1];

        const label charge = cloud_.constProps(id9).charge();

        if (charge != -1)
        {
            FatalErrorIn("dissociationIonisationExchange::setProperties()")
                << "Second product of the atom ionisation reaction must be "
                << "an electron: " << ionisationProductMoleculesQ[1]
                << nl
                << exit(FatalError);
        }
    }
}


bool dissociationIonisationExchange::tryReactMolecules(const label typeIdP,
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


void dissociationIonisationExchange::reaction
(
    dsmcParcel& p,
    dsmcParcel& q,
    const DynamicList<label>& candidateList,
    const List<DynamicList<label> >& candidateSubList,
    const label candidateP,
    const List<label>& whichSubCell
)
{}


void dissociationIonisationExchange::reaction
(
    dsmcParcel& p,
    dsmcParcel& q
)
{
    label typeIdP = p.typeId();
    label typeIdQ = q.typeId();

    // if particle p is the molecule and q is the atom
    if (typeIdP == reactantIds_[0] && typeIdQ == reactantIds_[1])
    {
        relax_ = true;

        scalar totalReactionProbability = 0.0;
        scalarList reactionProbabilities(5, 0.0);

        scalar EcPQ = 0.0;
        scalar TColl = 0.0;
        label idP = -1;
        label deltaDissoIP = 0;
        label imaxP = 0;
        label iaP = 0;
        bool dissocReaction = false;
        bool ionisationReactionP = false;
        bool ionisationReactionQ = false;
        bool exchangeReaction = false;
        bool chargeExchangeReaction = false;

        if (!chargedMolecule_)
        {
            // firstly calculate dissociation probability (0 or 1).
            EcPQ = reactionProperties(cloud_, p, q).translationalEnergy() + reactionProperties(cloud_, p, q).EVibP();

            imaxP = EcPQ/(physicoChemical::k.value()*reactionProperties(cloud_, p, q).thetaVP());

            idP = reactionProperties(cloud_, p, q).thetaDP()/reactionProperties(cloud_, p, q).thetaVP();

            deltaDissoIP = imaxP - idP;

            if (deltaDissoIP > 0)
            {
                totalReactionProbability += 1.0;
                reactionProbabilities[0] = 1.0;
            }

            // Now, ionisation of the molecule (P)

            scalar ionisationEnergy =
                            cloud_.constProps(typeIdP).ionisationTemperature()
                            *physicoChemical::k.value();

            // calculate if an ionisation of species P is possible
            EcPQ = reactionProperties(cloud_, p, q).translationalEnergy() + reactionProperties(cloud_, p, q).EEleP();

            if ((EcPQ - ionisationEnergy) > VSMALL)
            {
                // Ionisation can occur
//                 totalReactionProbability += 1.0;
//                 reactionProbabilities[1] = 1.0;
            }
        }

        if (!chargedAtom_)
        {
            // Now, ionisation of the atom (Q)

            scalar ionisationEnergy =
                            cloud_.constProps(typeIdQ).ionisationTemperature()
                            *physicoChemical::k.value();

            // calculate if an ionisation of species Q is possible
            EcPQ = reactionProperties(cloud_, p, q).translationalEnergy() + reactionProperties(cloud_, p, q).EEleQ();

            if ((EcPQ - ionisationEnergy) > VSMALL)
            {
                // Ionisation can occur
//                 totalReactionProbability += 1.0;
//                 reactionProbabilities[2] = 1.0;
            }
        }

        EcPQ = reactionProperties(cloud_, p, q).translationalEnergy() + reactionProperties(cloud_, p, q).EVibP();

        scalar P_exch = 0.0;

        // Next, calculate exchange probability.

        TColl = (reactionProperties(cloud_, p, q).translationalEnergy()/(physicoChemical::k.value()))
                                                    /(2.5 - reactionProperties(cloud_, p, q).omegaPQ());

        scalar heatOfReactionExchJoules = heatOfReactionExch_
                                                *physicoChemical::k.value();

        scalar aDash = aCoeff_*(pow(2.5 - reactionProperties(cloud_, p, q).omegaPQ(), bCoeff_)
                            *exp(lgamma(2.5 - reactionProperties(cloud_, p, q).omegaPQ()))
                            /exp(lgamma(2.5 - reactionProperties(cloud_, p, q).omegaPQ() + bCoeff_)));

        scalar activationEnergy = (aDash*pow((TColl/273.0), bCoeff_)
                                    * mag(heatOfReactionExchJoules));

        if (heatOfReactionExchJoules < 0.0)
        {
            activationEnergy -= heatOfReactionExchJoules;
        }

        if (EcPQ > activationEnergy)
        {
            scalar summation = 0.0; 

            if (activationEnergy <
             cloud_.constProps(typeIdP).thetaV()[0]*physicoChemical::k.value())
            {
                summation = 1.0;
            }
            else
            {
                iaP = EcPQ / (physicoChemical::k.value()*reactionProperties(cloud_, p, q).thetaVP());

                for (label i = 0; i <= iaP; i++)
                {
                    summation += pow(
                                (1.0 -((i*physicoChemical::k.value()
                                *reactionProperties(cloud_, p, q).thetaVP())
                                /EcPQ)),
                                (1.5 - reactionProperties(cloud_, p, q).omegaPQ()));
                }
            }

            P_exch = pow((1.0 - (activationEnergy/EcPQ)),
                                    (1.5 - reactionProperties(cloud_, p, q).omegaPQ()))/summation;

            totalReactionProbability += P_exch;
            reactionProbabilities[3] = P_exch;
        }

        if (chargeExchange_)
        {
            // calculate charge exchange probability

            label maxElectronicLevelP =
                    cloud_.constProps(typeIdP).numberOfElectronicLevels();

            scalar heatOfReactionChargeExchJoules =
                    heatOfReactionChargeExchange_*physicoChemical::k.value();

            scalar EcP = reactionProperties(cloud_, p, q).translationalEnergy() + reactionProperties(cloud_, p, q).EEleP();

            scalar aDash = aCoeffCharge_*(pow(2.5 - reactionProperties(cloud_, p, q).omegaPQ(), bCoeffCharge_)
                                *exp(lgamma(2.5 - reactionProperties(cloud_, p, q).omegaPQ()))
                                /exp(lgamma(2.5 - reactionProperties(cloud_, p, q).omegaPQ() + bCoeffCharge_)));

            TColl = (reactionProperties(cloud_, p, q).translationalEnergy()/(physicoChemical::k.value()))
                            /(2.5 - reactionProperties(cloud_, p, q).omegaPQ());

            scalar activationEnergy = (aDash*pow((TColl/273.0), bCoeffCharge_)
                                    * mag(heatOfReactionChargeExchJoules));

    //         scalar activationEnergy = activationEnergy_
           // + (aCoeff_*pow((5000/273.0), bCoeff_)
           // * mag(heatOfReactionExchJoules));

            if (heatOfReactionChargeExchJoules < 0.0)
            {
                activationEnergy -= heatOfReactionChargeExchJoules;
            }

            if (EcP > activationEnergy)
            {
                label keyElectronicLevel = -1;

                for (label i = 0; i < maxElectronicLevelP; i++)
                {
                    scalar electronicEnergy = reactionProperties(cloud_, p, q).EEListP()[i];

                    if (electronicEnergy > activationEnergy)
                    {
                        break;
                    }

                    keyElectronicLevel++;
                }

                EcP = reactionProperties(cloud_, p, q).translationalEnergy() + reactionProperties(cloud_, p, q).EEleP(); /*+heatOfReactionExchJoules*/

                label trialELevel = cloud_.postCollisionElectronicEnergyLevel
                                (
                                    EcP,
                                    maxElectronicLevelP,
                                    reactionProperties(cloud_, p, q).omegaPQ(),
                                    reactionProperties(cloud_, p, q).EEListP(),
                                    reactionProperties(cloud_, p, q).gListP()
                                );

                if (trialELevel == keyElectronicLevel)
                {
                    scalar prob = 0.0;

                    label nPossStates = 0;

                    if (maxElectronicLevelP == 1)
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

                    label nState = ceil(cloud_.rndGen().scalar01()
                                        *(nPossStates));
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
                        summation += reactionProperties(cloud_, p, q).gListP()[i]*pow((EcP - reactionProperties(cloud_, p, q).EEListP()[i]),
                                                        1.5 - reactionProperties(cloud_, p, q).omegaPQ());
                    }

                    prob = (reactionProperties(cloud_, p, q).gListP()[trialELevel]
                            *pow((EcP - reactionProperties(cloud_, p, q).EEListP()[trialELevel]), 1.5 - reactionProperties(cloud_, p, q).omegaPQ()))
                            /summation;

                    if (prob > cloud_.rndGen().scalar01())
                    {
                        // Charge exchange can occur
                        totalReactionProbability += prob;
                        reactionProbabilities[4] = prob;
                    }
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
                            // Dissociation is to occur
                            dissocReaction = true;
                            break;
                        }
                        if (i == 1)
                        {
                            // Molecule ionisation reaction is to occur
                            ionisationReactionP = true;
                            break;
                        }
                        if (i == 2)
                        {
                            // Atom ionisation reaction is to occur
                            ionisationReactionQ = true;
                            break;
                        }
                        if (i == 3)
                        {
                            // Exchange reaction is to occur
                            exchangeReaction = true;
                            break;
                        }
                        if (i == 4)
                        {
                            // Exchange reaction is to occur
                            chargeExchangeReaction = true;
                            break;
                        }
                    }
                }
            }
        }

        // Perform a dissociation reaction
        if (dissocReaction)
        {
            nTotDissociationReactions_++;
            nDissociationReactionsPerTimeStep_++;

            if (allowSplitting_)
            {
                relax_ = false;

                reactionProperties(cloud_, p, q).dissociateP
                (
                    heatOfReactionDiss_,
                    dissociationProductIds_,
                    p,
                    q
                );
            }
        }

        if (ionisationReactionP)
        {
            nTotIonisationReactionsP_++;
            nIonisationReactionsPPerTimeStep_++;

            if (allowSplitting_)
            {
                relax_ = false;
                
                reactionProperties(cloud_, p, q).ioniseP
                (
                    heatOfReactionIonP_,
                    ionisationProductsIdsP_,
                    p,
                    q
                );
            }
        }

        if (ionisationReactionQ)
        {
           nTotIonisationReactionsQ_++;
           nIonisationReactionsQPerTimeStep_++;

           if (allowSplitting_)
            {
                relax_ = false;
                
                reactionProperties(cloud_, p, q).ioniseQ
                (
                    heatOfReactionIonQ_,
                    ionisationProductsIdsQ_,
                    p,
                    q
                );

            }
        }

        // Perform exchange reaction
        if (exchangeReaction)
        {
            nTotExchangeReactions_++;
            nExchangeReactionsPerTimeStep_++;

            if (allowSplitting_)
            {
                relax_ = false;

                const label typeIdMol = exchangeProductIds_[0];
                const label typeIdAtom = exchangeProductIds_[1];

                // change species properties

                const scalar& mPExch = cloud_.constProps(typeIdAtom).mass();
                const scalar& mQExch = cloud_.constProps(typeIdMol).mass();

                scalar mRExch = mPExch*mQExch/(mPExch + mQExch);
                
                scalar translationalEnergy = reactionProperties(cloud_, p, q).translationalEnergy();

                translationalEnergy += (reactionProperties(cloud_, p, q).ERotP() + reactionProperties(cloud_, p, q).EVibP() + reactionProperties(cloud_, p, q).EEleP() 
                                        + reactionProperties(cloud_, p, q).EEleQ() + heatOfReactionExchJoules);

                scalar relVelExchMol = sqrt((2.0*translationalEnergy)/mRExch);

                // Variable Hard Sphere collision part for collision of molecules

                scalar cosTheta = 2.0*cloud_.rndGen().scalar01() - 1.0;

                scalar sinTheta = sqrt(1.0 - cosTheta*cosTheta);

                scalar phi = twoPi*cloud_.rndGen().scalar01();

                vector postCollisionRelU =
                    relVelExchMol
                    *vector
                        (
                            cosTheta,
                            sinTheta*cos(phi),
                            sinTheta*sin(phi)
                        );

                vector UP = reactionProperties(cloud_, p, q).Ucm() + (postCollisionRelU*mQExch/(mPExch + mQExch));
                vector UQ = reactionProperties(cloud_, p, q).Ucm() - (postCollisionRelU*mPExch/(mPExch + mQExch));

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
                labelList vibLevelAtom(0, 0);
                
                p.position() = positionP;
                p.U() = UP;
                p.RWF() = RWFp;
                p.ERot() = 0.0;
                p.ELevel() = 0;
                p.cell() = cell;
                p.tetFace() = tetFace;
                p.tetPt() = tetPt;
                p.typeId() = typeIdAtom;
                p.newParcel() = 0;
                p.vibLevel() = vibLevelAtom;
                
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
                labelList vibLevelMol(1, 0);
                
                q.position() = positionQ;
                q.U() = UQ;
                q.RWF() = RWFq;
                q.ERot() = 0.0;
                q.ELevel() = 0;
                q.cell() = cell;
                q.tetFace() = tetFace;
                q.tetPt() = tetPt;
                q.typeId() = typeIdMol;
                q.newParcel() = 0;
                q.vibLevel() = vibLevelMol;
            }
        }

        if (chargeExchangeReaction)
        {
            relax_ = false;

            nTotChargeExchangeReactions_++;
            nChargeExchangeReactionsPerTimeStep_++;

            if (allowSplitting_)
            {
                
                const label typeIdMol = chargeExchangeProductIds_[0];
                const label typeIdAtom = chargeExchangeProductIds_[1];
                
                scalar heatOfReactionChargeExchJoules =
                    heatOfReactionChargeExchange_*physicoChemical::k.value();

                scalar translationalEnergy = reactionProperties(cloud_, p, q).translationalEnergy();

                translationalEnergy += heatOfReactionChargeExchJoules;

                translationalEnergy += reactionProperties(cloud_, p, q).ERotP() + reactionProperties(cloud_, p, q).EVibP() + reactionProperties(cloud_, p, q).EEleP()
                                        + reactionProperties(cloud_, p, q).ERotQ() + reactionProperties(cloud_, p, q).EVibQ() + reactionProperties(cloud_, p, q).EEleQ();
                                        
                const scalar& mPExch = cloud_.constProps(typeIdAtom).mass();
                const scalar& mQExch = cloud_.constProps(typeIdMol).mass();

                scalar mRExch = mPExch*mQExch/(mPExch + mQExch);

                scalar relVel = sqrt((2.0*translationalEnergy)/mRExch);

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

                vector UP = reactionProperties(cloud_, p, q).Ucm() + (postCollisionRelU*mQExch/(mPExch + mQExch));
                vector UQ = reactionProperties(cloud_, p, q).Ucm() - (postCollisionRelU*mPExch/(mPExch + mQExch));
                
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
                labelList vibLevelAtom(0, 0);
                
                p.position() = positionP;
                p.U() = UP;
                p.RWF() = RWFp;
                p.ERot() = 0.0;
                p.ELevel() = 0;
                p.cell() = cell;
                p.tetFace() = tetFace;
                p.tetPt() = tetPt;
                p.typeId() = typeIdAtom;
                p.newParcel() = 0;
                p.vibLevel() = vibLevelAtom;

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
                labelList vibLevelMol(1, 0);
                
                q.position() = positionQ;
                q.U() = UQ;
                q.RWF() = RWFq;
                q.ERot() = 0.0;
                q.ELevel() = 0;
                q.cell() = cell;
                q.tetFace() = tetFace;
                q.tetPt() = tetPt;
                q.typeId() = typeIdMol;
                q.newParcel() = 0;
                q.vibLevel() = vibLevelMol;
            }
        }
    }

    // if q is the molecule and p is the atom...
    if (typeIdP == reactantIds_[1] && typeIdQ == reactantIds_[0])
    {
        relax_ = true;

        scalar totalReactionProbability = 0.0;
        scalarList reactionProbabilities(5, 0.0);

        scalar EcPQ = 0.0;
        scalar TColl = 0.0;
        label idQ = -1;
        label deltaDissoIQ = 0;
        label imaxQ = 0;
        label iaQ = 0;
        bool dissocReaction = false;
        bool ionisationReactionP = false;
        bool ionisationReactionQ = false;
        bool exchangeReaction = false;
        bool chargeExchangeReaction = false;

        if (!chargedMolecule_)
        {
            // firstly calculate dissociation probability (0 or 1).
            EcPQ = reactionProperties(cloud_, p, q).translationalEnergy() + reactionProperties(cloud_, p, q).EVibQ();

            imaxQ = EcPQ/(physicoChemical::k.value()*reactionProperties(cloud_, p, q).thetaVQ());

            idQ = reactionProperties(cloud_, p, q).thetaDQ()/reactionProperties(cloud_, p, q).thetaVQ();

            deltaDissoIQ = imaxQ - idQ;

            if (deltaDissoIQ > 0)
            {
                totalReactionProbability += 1.0;
                reactionProbabilities[0] = 1.0;
            }

            scalar ionisationEnergy =
                        cloud_.constProps(typeIdQ).ionisationTemperature()
                        *physicoChemical::k.value();

            // calculate if an ionisation of species Q is possible
            EcPQ = reactionProperties(cloud_, p, q).translationalEnergy() + reactionProperties(cloud_, p, q).EEleQ();

            if ((EcPQ - ionisationEnergy) > VSMALL)
            {
                // Ionisation can occur
//                 totalReactionProbability += 1.0;
//                 reactionProbabilities[1] = 1.0;
            }
        }

        if (!chargedAtom_)
        {
            scalar ionisationEnergy =
                        cloud_.constProps(typeIdP).ionisationTemperature()
                        *physicoChemical::k.value();

            EcPQ = reactionProperties(cloud_, p, q).translationalEnergy() + reactionProperties(cloud_, p, q).EEleP();

            if ((EcPQ - ionisationEnergy) > VSMALL)
            {
                // Ionisation can occur
//                 totalReactionProbability += 1.0;
//                 reactionProbabilities[2] = 1.0;
            }
        }

        EcPQ = reactionProperties(cloud_, p, q).translationalEnergy() + reactionProperties(cloud_, p, q).EVibQ();

        scalar P_exch = 0.0;

        TColl = (reactionProperties(cloud_, p, q).translationalEnergy()/(physicoChemical::k.value()))
                            /(2.5 - reactionProperties(cloud_, p, q).omegaPQ());

        scalar heatOfReactionExchJoules = heatOfReactionExch_
                                            *physicoChemical::k.value();

        scalar aDash = aCoeff_*(pow(2.5 - reactionProperties(cloud_, p, q).omegaPQ(), bCoeff_)
                        *exp(lgamma(2.5 - reactionProperties(cloud_, p, q).omegaPQ()))
                        /exp(lgamma(2.5 - reactionProperties(cloud_, p, q).omegaPQ() + bCoeff_)));

        scalar activationEnergy = (aDash*pow((TColl/273.0), bCoeff_)
                                    * mag(heatOfReactionExchJoules));

        if (heatOfReactionExchJoules < 0.0)
        {
            activationEnergy -= heatOfReactionExchJoules;
        }

        if (EcPQ > activationEnergy)
        {
            scalar summation = 0.0;

            if (activationEnergy <
             cloud_.constProps(typeIdQ).thetaV()[0]*physicoChemical::k.value())
            {
                // this refers to the first sentence in Bird's QK paper after
                // Eq.(12)
                summation = 1.0;
            }
            else
            {
                iaQ = EcPQ / (physicoChemical::k.value()*reactionProperties(cloud_, p, q).thetaVQ());

                for (label i = 0; i <= iaQ; i++)
                {
                    summation += pow((1.0 -
                                    ((i*physicoChemical::k.value()
                                    *reactionProperties(cloud_, p, q).thetaVQ())
                                    /EcPQ)),
                                    (1.5 - reactionProperties(cloud_, p, q).omegaPQ()));
                }
            }

            P_exch = pow((1.0 - (activationEnergy/EcPQ)),
                         (1.5 - reactionProperties(cloud_, p, q).omegaPQ()))/summation;

            totalReactionProbability += P_exch;
            reactionProbabilities[3] = P_exch;
        }

        if (chargeExchange_)
        {
            label maxElectronicLevelQ =
                    cloud_.constProps(typeIdQ).numberOfElectronicLevels();

            scalar heatOfReactionChargeExchJoules =
                    heatOfReactionChargeExchange_*physicoChemical::k.value();

            scalar EcQ = reactionProperties(cloud_, p, q).translationalEnergy() + reactionProperties(cloud_, p, q).EEleQ();

            scalar aDash = aCoeffCharge_*(pow(2.5 - reactionProperties(cloud_, p, q).omegaPQ(), bCoeffCharge_)
                                *exp(lgamma(2.5 - reactionProperties(cloud_, p, q).omegaPQ()))
                                /exp(lgamma(2.5 - reactionProperties(cloud_, p, q).omegaPQ() +
                                bCoeffCharge_)));

            TColl = (reactionProperties(cloud_, p, q).translationalEnergy()/(physicoChemical::k.value()))
                        /(2.5 - reactionProperties(cloud_, p, q).omegaPQ());

            scalar activationEnergy = (aDash*pow((TColl/273.0), bCoeffCharge_)
                                * mag(heatOfReactionChargeExchJoules));

    //         scalar activationEnergy = (aCoeff_*pow((5000/273.0), bCoeff_)
                            //* mag(heatOfReactionExchJoules));

            if (heatOfReactionChargeExchJoules < 0.0)
            {
                activationEnergy -= heatOfReactionChargeExchJoules;
            }

            if (EcQ > activationEnergy)
            {
                label keyElectronicLevel = -1;

                for (label i = 0; i < maxElectronicLevelQ; i++)
                {
                    scalar electronicEnergy = reactionProperties(cloud_, p, q).EEListQ()[i];

                    if (electronicEnergy > activationEnergy)
                    {
                        break;
                    }

                    keyElectronicLevel++;
                }

                EcQ = reactionProperties(cloud_, p, q).translationalEnergy() + reactionProperties(cloud_, p, q).EEleQ(); /*+heatOfReactionExchJoules*/

                label trialELevel = cloud_.postCollisionElectronicEnergyLevel
                                (
                                    EcQ,
                                    maxElectronicLevelQ,
                                    reactionProperties(cloud_, p, q).omegaPQ(),
                                    reactionProperties(cloud_, p, q).EEListQ(),
                                    reactionProperties(cloud_, p, q).gListQ()
                                );

                if (trialELevel == keyElectronicLevel)
                {
                    scalar prob = 0.0;

                    label nPossStates = 0;

                    if (maxElectronicLevelQ == 1)
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

                    label nState = ceil(cloud_.rndGen().scalar01()
                                            *(nPossStates));
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

                    // Calculate the probability of it occuring

                    scalar summation = 0.0;

                    for (label i = 0; i <= nLevel; i++)
                    {
                        summation += reactionProperties(cloud_, p, q).gListQ()[i]
                                    *pow((EcQ - reactionProperties(cloud_, p, q).EEListQ()[i]), 1.5 - reactionProperties(cloud_, p, q).omegaPQ());
                    }

                    prob = (reactionProperties(cloud_, p, q).gListQ()[trialELevel]
                            *pow((EcQ - reactionProperties(cloud_, p, q).EEListQ()[trialELevel]),
                                1.5 - reactionProperties(cloud_, p, q).omegaPQ()))/summation;

                    if (prob > cloud_.rndGen().scalar01())
                    {
                        // Charge exchange can occur
                        totalReactionProbability += prob;
                        reactionProbabilities[4] = prob;
                    }
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
                            // Dissociation is to occur
                            dissocReaction = true;
                            break;
                        }
                        if (i == 1)
                        {
                            // Molecule ionisation reaction is to occur
                            ionisationReactionP = true;
                            break;
                        }
                        if (i == 2)
                        {
                            // Atom ionisation reaction is to occur
                            ionisationReactionQ = true;
                            break;
                        }
                        if (i == 3)
                        {
                            // Exchange reaction is to occur
                            exchangeReaction = true;
                            break;
                        }
                        if (i == 4)
                        {
                            // Exchange reaction is to occur
                            chargeExchangeReaction = true;
                            break;
                        }
                    }
                }
            }
        }

        // Perform dissociation reaction
        if (dissocReaction)
        {
            nTotDissociationReactions_++;
            nDissociationReactionsPerTimeStep_++;

            if (allowSplitting_)
            {
                relax_ = false;
                
                reactionProperties(cloud_, p, q).dissociateQ
                (
                    heatOfReactionDiss_,
                    dissociationProductIds_,
                    p,
                    q
                );
            }
        }

        if (ionisationReactionP)
        {
            // Molecule ionisation (Q is the molecule, P is used
            // for measurement purposes)
            nTotIonisationReactionsP_++;
            nIonisationReactionsPPerTimeStep_++;

            if (allowSplitting_)
            {
                relax_ = false;
                
                reactionProperties(cloud_, p, q).ioniseQ
                (
                    heatOfReactionIonP_,
                    ionisationProductsIdsP_,
                    p,
                    q
                );
            }
        }

        if (ionisationReactionQ)
        {
           // P is the atom (Q is used for measurement purposes)
           nTotIonisationReactionsQ_++;
           nIonisationReactionsQPerTimeStep_++;

           if (allowSplitting_)
           {
                relax_ = false;
                
                reactionProperties(cloud_, p, q).ioniseP
                (
                    heatOfReactionIonQ_,
                    ionisationProductsIdsQ_,
                    p,
                    q
                );
            }
        }

        // Perform exchange reaction
        if (exchangeReaction)
        {
            nTotExchangeReactions_++;
            nExchangeReactionsPerTimeStep_++;

            if (allowSplitting_)
            {
                relax_ = false;

                const label typeIdMol = exchangeProductIds_[0];
                const label typeIdAtom = exchangeProductIds_[1];

                // change species properties

                const scalar& mQExch = cloud_.constProps(typeIdAtom).mass();
                const scalar& mPExch = cloud_.constProps(typeIdMol).mass();

                scalar mRExch = mPExch*mQExch/(mPExch + mQExch);
                
                scalar translationalEnergy = reactionProperties(cloud_, p, q).translationalEnergy();

                translationalEnergy += reactionProperties(cloud_, p, q).ERotQ() + reactionProperties(cloud_, p, q).EVibQ()
                                        + reactionProperties(cloud_, p, q).EEleQ() + reactionProperties(cloud_, p, q).EEleP()
                                        + heatOfReactionExchJoules;

                scalar relVelExchMol = sqrt((2.0*translationalEnergy)/mRExch);

                // Variable Hard Sphere collision part for collision of
                // molecules

                scalar cosTheta = 2.0*cloud_.rndGen().scalar01() - 1.0;

                scalar sinTheta = sqrt(1.0 - cosTheta*cosTheta);

                scalar phi = twoPi*cloud_.rndGen().scalar01();

                vector postCollisionRelU =
                    relVelExchMol
                    *vector
                        (
                            cosTheta,
                            sinTheta*cos(phi),
                            sinTheta*sin(phi)
                        );

                vector UP = reactionProperties(cloud_, p, q).Ucm() + (postCollisionRelU*mQExch/(mPExch + mQExch));
                vector UQ = reactionProperties(cloud_, p, q).Ucm() - (postCollisionRelU*mPExch/(mPExch + mQExch));

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
                labelList vibLevelMol(1, 0);
                
                p.position() = positionP;
                p.U() = UP;
                p.RWF() = RWFp;
                p.ERot() = 0.0;
                p.ELevel() = 0;
                p.cell() = cell;
                p.tetFace() = tetFace;
                p.tetPt() = tetPt;
                p.typeId() = typeIdMol;
                p.newParcel() = 0;
                p.vibLevel() = vibLevelMol;
                
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
                labelList vibLevelAtom(0, 0);
                
                q.position() = positionQ;
                q.U() = UQ;
                q.RWF() = RWFq;
                q.ERot() = 0.0;
                q.ELevel() = 0;
                q.cell() = cell;
                q.tetFace() = tetFace;
                q.tetPt() = tetPt;
                q.typeId() = typeIdAtom;
                q.newParcel() = 0;
                q.vibLevel() = vibLevelAtom;
            }
        }

        if (chargeExchangeReaction)
        {
            nTotChargeExchangeReactions_++;
            nChargeExchangeReactionsPerTimeStep_++;

            if (allowSplitting_)
            {
                relax_ = false;
                
                const label typeIdMol = chargeExchangeProductIds_[0];
                const label typeIdAtom = chargeExchangeProductIds_[1];
                
                scalar heatOfReactionChargeExchJoules =
                    heatOfReactionChargeExchange_*physicoChemical::k.value();

                scalar translationalEnergy = reactionProperties(cloud_, p, q).translationalEnergy();

                translationalEnergy += heatOfReactionChargeExchJoules;

                translationalEnergy += reactionProperties(cloud_, p, q).ERotP() + reactionProperties(cloud_, p, q).EVibP() + reactionProperties(cloud_, p, q).EEleP()
                                        + reactionProperties(cloud_, p, q).ERotQ() + reactionProperties(cloud_, p, q).EVibQ() + reactionProperties(cloud_, p, q).EEleQ();
                                        
                const scalar& mPExch = cloud_.constProps(typeIdMol).mass();
                const scalar& mQExch = cloud_.constProps(typeIdAtom).mass();

                scalar mRExch = mPExch*mQExch/(mPExch + mQExch);

                scalar relVel = sqrt((2.0*translationalEnergy)/mRExch);

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

                vector UP = reactionProperties(cloud_, p, q).Ucm() + (postCollisionRelU*mQExch/(mPExch + mQExch));
                vector UQ = reactionProperties(cloud_, p, q).Ucm() - (postCollisionRelU*mPExch/(mPExch + mQExch));
                
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
                labelList vibLevelMol(1, 0);
                
                p.position() = positionP;
                p.U() = UP;
                p.RWF() = RWFp;
                p.ERot() = 0.0;
                p.ELevel() = 0;
                p.cell() = cell;
                p.tetFace() = tetFace;
                p.tetPt() = tetPt;
                p.typeId() = typeIdMol;
                p.newParcel() = 0;
                p.vibLevel() = vibLevelMol;

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
                labelList vibLevelAtom(0, 0);
                
                q.position() = positionQ;
                q.U() = UQ;
                q.RWF() = RWFq;
                q.ERot() = 0.0;
                q.ELevel() = 0;
                q.cell() = cell;
                q.tetFace() = tetFace;
                q.tetPt() = tetPt;
                q.typeId() = typeIdAtom;
                q.newParcel() = 0;
                q.vibLevel() = vibLevelAtom;
            }
        }
    }
}


void dissociationIonisationExchange::reactExchangeMolecule
(
    dsmcParcel& p,
    label newTypeId,
    const label newEVibLevel,
    const scalar newERot,
    const vector& newU
)
{
    p.vibLevel() = newEVibLevel;
    p.ERot() = newERot;
    p.U() = newU;
    p.typeId() = newTypeId;
}


void dissociationIonisationExchange::reactExchangeAtom
(
    dsmcParcel& p,
    label newTypeId,
    const vector& newU
)
{
    p.U() = newU;
    p.typeId() = newTypeId;
}


void  dissociationIonisationExchange::outputResults(const label counterIndex)
{
    if (writeRatesToTerminal_ == true)
    {
        const List<DynamicList<dsmcParcel*> >& cellOccupancy
            = cloud_.cellOccupancy();

        volume_ = 0.0;

        forAll(cellOccupancy, c)
        {
            volume_ += mesh_.cellVolumes()[c];
        }

        List<label> mols (2, 0);
        scalar volume = volume_;
        label nTotExchangeReactions = nTotExchangeReactions_;
        label nTotChargeExchangeReactions = nTotChargeExchangeReactions_;
        label nTotDissociationReactions = nTotDissociationReactions_;
        label nTotIonisationReactionsP = nTotIonisationReactionsP_;
        label nTotIonisationReactionsQ = nTotIonisationReactionsQ_;

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
            reduce(nTotExchangeReactions, sumOp<label>());
            reduce(nTotChargeExchangeReactions, sumOp<label>());
            reduce(nTotDissociationReactions, sumOp<label>());
            reduce(nTotIonisationReactionsP, sumOp<label>());
            reduce(nTotIonisationReactionsQ, sumOp<label>());
        }

        numberDensities_[0] = (mols[0]*cloud().nParticle())/volume;
        numberDensities_[1] = (mols[1]*cloud().nParticle())/volume;

        const scalar deltaT = mesh_.time().deltaT().value();

        word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
        word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];

        word exchangeProductMolA = cloud_.typeIdList()[exchangeProductIds_[0]];
        word exchangeProductMolB = cloud_.typeIdList()[exchangeProductIds_[1]];

        word dissociationProductMolA;
        word dissociationProductMolB;

        word moleculeIonisationProductMolA;
        word moleculeIonisationProductMolB;

        word atomIonisationProductMolA;
        word atomIonisationProductMolB;

        word chargeExchangeProductMolA;
        word chargeExchangeProductMolB;

        if (!chargedMolecule_)
        {
            dissociationProductMolA =
                    cloud_.typeIdList()[dissociationProductIds_[0]];
            dissociationProductMolB =
                    cloud_.typeIdList()[dissociationProductIds_[1]];

            moleculeIonisationProductMolA =
                    cloud_.typeIdList()[ionisationProductsIdsP_[0]];
            moleculeIonisationProductMolB =
                    cloud_.typeIdList()[ionisationProductsIdsP_[1]];
        }

        if (!chargedAtom_)
        {
            atomIonisationProductMolA =
                    cloud_.typeIdList()[ionisationProductsIdsQ_[0]];
            atomIonisationProductMolB =
                    cloud_.typeIdList()[ionisationProductsIdsQ_[1]];
        }

        if (chargeExchange_)
        {
            chargeExchangeProductMolA =
                    cloud_.typeIdList()[chargeExchangeProductIds_[0]];
            chargeExchangeProductMolB =
                    cloud_.typeIdList()[chargeExchangeProductIds_[1]];
        }

        if ((numberDensities_[0] > 0.0) && (numberDensities_[1] > 0.0))
        {
            scalar reactionRate1 = 0.0;
            scalar reactionRate2 = 0.0;
            scalar reactionRate3 = 0.0;
            scalar reactionRate4 = 0.0;
            scalar reactionRate5 = 0.0;

            reactionRate1 =
            (
                nTotExchangeReactions
                * cloud_.nParticle()
)/(counterIndex*deltaT*numberDensities_[0]
            * numberDensities_[1]*volume);

            Info << "Exchange reaction "
                    <<  reactantMolA << " + " << reactantMolB << " --> "
                    << exchangeProductMolA << " + " << exchangeProductMolB
                    << ", reaction rate = " << reactionRate1
                    << endl;


            if (!chargedMolecule_)
            {
                reactionRate2 =
                (
                    nTotDissociationReactions
                    * cloud_.nParticle()
)/(counterIndex*deltaT*numberDensities_[0]
                    * numberDensities_[1]*volume);

                Info << "Type II dissociation reaction "
                    <<  reactantMolA << " + " << reactantMolB << " --> "
                    << dissociationProductMolA << " + "
                    << dissociationProductMolB
                    <<  " + " << reactantMolB
                    << ", reaction rate = " << reactionRate2
                    << endl;


                reactionRate3 =
                (
                    nTotIonisationReactionsP
                    * cloud_.nParticle()
)/(counterIndex*deltaT*numberDensities_[0]
                    * numberDensities_[1]*volume);

                Info << "Ionisation reaction "
                    <<  reactantMolA << " + " << reactantMolB << " --> "
                    << moleculeIonisationProductMolA << " + "
                    << moleculeIonisationProductMolB
                    <<  " + " << reactantMolB
                    << ", reaction rate = " << reactionRate3
                    << endl;
            }

            if (!chargedAtom_)
            {
                reactionRate4 =
                (
                    nTotIonisationReactionsQ
                    * cloud_.nParticle()
)/(counterIndex*deltaT*numberDensities_[0]
                    * numberDensities_[1]*volume);

                Info << "Ionisation reaction "
                    <<  reactantMolA << " + " << reactantMolB << " --> "
                    <<  reactantMolA << " + "
                    << atomIonisationProductMolA << " + "
                    << atomIonisationProductMolB
                    << ", reaction rate = " << reactionRate4
                    << endl;

            }

            if (chargeExchange_)
            {
                reactionRate5 =
                (
                    nTotChargeExchangeReactions
                    * cloud_.nParticle()
)/(counterIndex*deltaT*numberDensities_[0]
                    * numberDensities_[1]*volume);

               Info << "Charge exchange reaction "
                    <<  reactantMolA << " + " << reactantMolB << " --> "
                    << chargeExchangeProductMolA << " + "
                    << chargeExchangeProductMolB
                    << ", reaction rate = " << reactionRate5
                    << endl;
            }
        }
    }
    else
    {
        word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
        word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];

        word exchangeProductMolA = cloud_.typeIdList()[exchangeProductIds_[0]];
        word exchangeProductMolB = cloud_.typeIdList()[exchangeProductIds_[1]];

        word dissociationProductMolA;
        word dissociationProductMolB;

        word moleculeIonisationProductMolA;
        word moleculeIonisationProductMolB;

        if (!chargedMolecule_)
        {
            dissociationProductMolA =
                    cloud_.typeIdList()[dissociationProductIds_[0]];
            dissociationProductMolB =
                    cloud_.typeIdList()[dissociationProductIds_[1]];

            moleculeIonisationProductMolA =
                    cloud_.typeIdList()[ionisationProductsIdsP_[0]];
            moleculeIonisationProductMolB =
                    cloud_.typeIdList()[ionisationProductsIdsP_[1]];
        }

        word atomIonisationProductMolA;
        word atomIonisationProductMolB;

        if (!chargedAtom_)
        {
            atomIonisationProductMolA =
                    cloud_.typeIdList()[ionisationProductsIdsQ_[0]];
            atomIonisationProductMolB =
                    cloud_.typeIdList()[ionisationProductsIdsQ_[1]];
        }

        word chargeExchangeProductMolA;
        word chargeExchangeProductMolB;

        if (chargeExchange_)
        {
            chargeExchangeProductMolA =
                    cloud_.typeIdList()[chargeExchangeProductIds_[0]];
            chargeExchangeProductMolB =
                    cloud_.typeIdList()[chargeExchangeProductIds_[1]];
        }

        label nTotExchangeReactions = nTotExchangeReactions_;
        label nTotChargeExchangeReactions = nTotChargeExchangeReactions_;
        label nTotDissociationReactions = nTotDissociationReactions_;
        label nTotIonisationReactionsP = nTotIonisationReactionsP_;
        label nTotIonisationReactionsQ = nTotIonisationReactionsQ_;

        label nDissociationReactionsPerTimeStep =
                nDissociationReactionsPerTimeStep_;
        label nIonisationReactionsPPerTimeStep =
                nIonisationReactionsPPerTimeStep_;
        label nIonisationReactionsQPerTimeStep =
                nIonisationReactionsQPerTimeStep_;
        label nExchangeReactionsPerTimeStep =
                nExchangeReactionsPerTimeStep_;
        label nChargeExchangeReactionsPerTimeStep =
                nChargeExchangeReactionsPerTimeStep_;

        // Parallel communication
        if (Pstream::parRun())
        {
            reduce(nTotExchangeReactions, sumOp<label>());
            reduce(nTotChargeExchangeReactions, sumOp<label>());
            reduce(nTotDissociationReactions, sumOp<label>());
            reduce(nTotIonisationReactionsP, sumOp<label>());
            reduce(nTotIonisationReactionsQ, sumOp<label>());

            reduce(nDissociationReactionsPerTimeStep, sumOp<label>());
            reduce(nIonisationReactionsPPerTimeStep, sumOp<label>());
            reduce(nIonisationReactionsQPerTimeStep, sumOp<label>());
            reduce(nExchangeReactionsPerTimeStep, sumOp<label>());
            reduce(nChargeExchangeReactionsPerTimeStep, sumOp<label>());
        }

        if (nTotExchangeReactions > VSMALL)
        {
            Info << "Exchange reaction "
                    <<  reactantMolA << " + " << reactantMolB << " --> "
                    << exchangeProductMolA << " + " << exchangeProductMolB
                    << " is active, nReactions this time step = "
                    << nExchangeReactionsPerTimeStep << endl;
        }

        if (nTotDissociationReactions > VSMALL)
        {
            Info  << "Type II dissociation reaction "
                <<  reactantMolA << " + " << reactantMolB << " --> "
                << dissociationProductMolA << " + " << dissociationProductMolB
                <<  " + " << reactantMolB
                << " is active, nReactions this time step = "
                << nDissociationReactionsPerTimeStep << endl;
        }

        if (nTotIonisationReactionsP > VSMALL)
        {
            Info  << "Ionisation reaction "
                <<  reactantMolA << " + " << reactantMolB << " --> "
                << moleculeIonisationProductMolA << " + "
                << moleculeIonisationProductMolB
                <<  " + " << reactantMolB
                << " is active, nReactions this time step = "
                << nIonisationReactionsPPerTimeStep << endl;
        }

        if (nTotIonisationReactionsQ > VSMALL)
        {
            Info  << "Ionisation reaction "
                <<  reactantMolA << " + " << reactantMolB << " --> "
                <<  reactantMolA << " + "
                << atomIonisationProductMolA << " + "
                << atomIonisationProductMolB
                << " is active, nReactions this time step = "
                << nIonisationReactionsQPerTimeStep << endl;
        }

        if (nTotChargeExchangeReactions > VSMALL)
        {
            Info  << "Charge exchange reaction "
                <<  reactantMolA << " + " << reactantMolB << " --> "
                << chargeExchangeProductMolA << " + "
                << chargeExchangeProductMolB
                << " is active, nReactions this time step = "
                << nChargeExchangeReactionsPerTimeStep << endl;
        }
    }

    nDissociationReactionsPerTimeStep_ = 0.0;
    nIonisationReactionsPPerTimeStep_ = 0.0;
    nIonisationReactionsQPerTimeStep_ = 0.0;
    nExchangeReactionsPerTimeStep_ = 0.0;
    nChargeExchangeReactionsPerTimeStep_ = 0.0;
}


bool dissociationIonisationExchange::relax() const
{
    return relax_;
}

} // End namespace Foam

// ************************************************************************* //
