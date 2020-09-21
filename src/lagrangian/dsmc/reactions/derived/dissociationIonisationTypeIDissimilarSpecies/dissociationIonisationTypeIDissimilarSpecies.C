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

#include "dissociationIonisationTypeIDissimilarSpecies.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(dissociationIonisationTypeIDissimilarSpecies, 0);

addToRunTimeSelectionTable
(
    dsmcReaction,
    dissociationIonisationTypeIDissimilarSpecies,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

dissociationIonisationTypeIDissimilarSpecies::
dissociationIonisationTypeIDissimilarSpecies
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
    heatOfReactionDissociationCD_(),
    heatOfReactionIonisationCD_(),
    nTotABDissociationReactions_(0),
    nABDissociationReactionsPerTimeStep_(0),
    nTotABIonisationReactions_(0),
    nABIonisationReactionsPerTimeStep_(0),
    nTotCDDissociationReactions_(0),
    nCDDissociationReactionsPerTimeStep_(0),
    nTotCDIonisationReactions_(0),
    nCDIonisationReactionsPerTimeStep_(0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dissociationIonisationTypeIDissimilarSpecies::
~dissociationIonisationTypeIDissimilarSpecies()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dissociationIonisationTypeIDissimilarSpecies::initialConfiguration()
{
    setCommonReactionProperties();
    setProperties();
}


void dissociationIonisationTypeIDissimilarSpecies::setProperties()
{
    // check that reactants are 'MOLECULES' (not 'ATOMS')

    if (rDof1_ < 1)
    {
        FatalErrorIn("dissociationIonisationTypeIDissimilarSpecies::setProperties()")
            << "Reactant must be a molecule (not an atom): "
            << reactants_[0]
            << nl
            << exit(FatalError);
    }

    if (vDof1_ > 1)
    {
        FatalErrorIn("dissociationIonisationTypeIDissimilarSpecies::setProperties()")
            << "Reactions are currently only implemented for "
            << "monatomic and diatomic species"
            << " This is a polyatomic:" << reactants_[0]
            << nl
            << exit(FatalError);
    }
    
    if (rDof2_ < 1)
    {
        FatalErrorIn("dissociationIonisationTypeIDissimilarSpecies::setProperties()")
            << "Reactant must be a molecule (not an atom): "
            << reactants_[1]
            << nl
            << exit(FatalError);
    }

    if (vDof2_ > 1)
    {
        FatalErrorIn("dissociationIonisationTypeIDissimilarSpecies::setProperties()")
            << "Reactions are currently only implemented for "
            << "monatomic and diatomic species"
            << " This is a polyatomic:" << reactants_[1]
            << nl
            << exit(FatalError);
    }

    // reading in dissociation products

    List<List<word> > dissociationProducts
        (propsDict_.lookup("productsOfDissociatedMolecule"));

    if (dissociationProducts.size() != reactantIds_.size())
    {
        FatalErrorIn("dissociationIonisationTypeIDissimilarSpecies::setProperties()")
            << "number of reactant molecules to be dissociated = "
            << reactantIds_.size()
            << " is not the same as the number of products = "
            << dissociationProducts.size()
            << exit(FatalError);
    }


    dissociationProducts_.setSize(dissociationProducts.size());

    forAll(dissociationProducts_, r)
    {
        const List<word>& productsForDiss = dissociationProducts[r];

        if (productsForDiss.size() != 2)
        {
            FatalErrorIn("dissociationIonisationTypeIDissimilarSpecies::setProperties()")
                << "There should be two products "
                << "(for the dissociating molecule "
                << reactants_[r] << "), instead of "
                << productsForDiss.size() << nl
                << exit(FatalError);
        }

        dissociationProducts_[r].setSize(productsForDiss.size(), -1);

        forAll(dissociationProducts_[r], p)
        {
            dissociationProducts_[r][p] = findIndex(cloud_.typeIdList(),
                                                    productsForDiss[p]);

            if (dissociationProducts_[r][p] == -1)
            {
                FatalErrorIn("dissociationIonisationTypeIDissimilarSpecies::setProperties()")
                    << "Cannot find type id: " << productsForDiss[p] << nl
                    << exit(FatalError);
            }
            
            // check that products are 'ATOMS' (not 'MOLECULES')
            label iD = dissociationProducts_[r][p];

            scalar rDof3 =
                    cloud_.constProps(iD).rotationalDegreesOfFreedom();

            if (rDof3 > 1)
            {
                FatalErrorIn("dissociationIonisationTypeIDissimilarSpecies::setProperties()")
                    << "Dissociation product must be an atom (not a molecule): "
                    << productsForDiss[r][p]
                    << nl
                    << exit(FatalError);
            }
        }
    }

    // reading in ionisation products

    List<List<word> > ionisationProducts
                (propsDict_.lookup("productsOfIonisedMolecule"));

    if (ionisationProducts.size() != reactantIds_.size())
    {
        FatalErrorIn("dissociationIonisationTypeIDissimilarSpecies::setProperties()")
            << "number of reactant molecules to be ionised = "
            << reactantIds_.size()
            << " is not the same as the number of products = "
            << ionisationProducts.size()
            << exit(FatalError);
    }


    ionisationProducts_.setSize(ionisationProducts.size());

    forAll(ionisationProducts_, r)
    {
        const List<word>& productsForIon = ionisationProducts[r];

        if (productsForIon.size() != 2)
        {
            FatalErrorIn("dissociationIonisationTypeIDissimilarSpecies::setProperties()")
                << "There should be two products (for the ionising molecule "
                << reactants_[r] << "), instead of "
                << productsForIon.size() << nl
                << exit(FatalError);
        }

        ionisationProducts_[r].setSize(productsForIon.size(), -1);

        forAll(ionisationProducts_[r], p)
        {
            ionisationProducts_[r][p] = findIndex(cloud_.typeIdList(),
                                                  productsForIon[p]);

            if (ionisationProducts_[r][p] == -1)
            {
                FatalErrorIn("dissociationIonisationTypeIDissimilarSpecies::setProperties()")
                    << "Cannot find type id: " << productsForIon[p] << nl
                    << exit(FatalError);
            }
        }

        // check ionisation product 2 is an electron
        forAll(ionisationProducts_[r], p)
        {
            if (p == 1)
            {
                const label charge =
                    cloud_.constProps(ionisationProducts_[r][p]).charge();

                if (charge != -1)
                {
                    FatalErrorIn("dissociationIonisationTypeIDissimilarSpecies::setProperties()")
                        << "Second ionisation product must be an electron: "
                        << productsForIon[p] << nl
                        << exit(FatalError);
                }
            }
        }
    }
    
    heatOfReactionDissociationAB_ =
        readScalar(propsDict_.lookup("heatOfReactionDissociationAB"));
    heatOfReactionIonisationAB_ =
        readScalar(propsDict_.lookup("heatOfReactionIonisationAB"));
    heatOfReactionDissociationCD_ =
        readScalar(propsDict_.lookup("heatOfReactionDissociationCD"));
    heatOfReactionIonisationCD_ =
        readScalar(propsDict_.lookup("heatOfReactionIonisationCD"));
}


bool dissociationIonisationTypeIDissimilarSpecies::
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


void dissociationIonisationTypeIDissimilarSpecies::reaction
(
    dsmcParcel& p,
    dsmcParcel& q,
    const DynamicList<label>& candidateList,
    const List<DynamicList<label> >& candidateSubList,
    const label candidateP,
    const List<label>& whichSubCell
)
{}


void dissociationIonisationTypeIDissimilarSpecies::reaction
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
        scalarList reactionProbabilities(4, 0.0);
        
        relax_ = true;

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
        label idP = reactionProperties(cloud_, p, q).charDissLevelP();
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

        scalar EcQP = 0.0;
        label idQ = reactionProperties(cloud_, p, q).charDissLevelQ();
        label imaxQ = 0;

        // calculate if a dissociation of species Q is possible
        EcQP = reactionProperties(cloud_, p, q).translationalEnergy() + reactionProperties(cloud_, p, q).EVibQ();

        imaxQ = EcQP/(physicoChemical::k.value()*reactionProperties(cloud_, p, q).thetaVQ());

        if (imaxQ - idQ > 0)
        {
            // Dissociation can occur
            totalReactionProbability += 1.0;
            reactionProbabilities[2] = 1.0;
        }

        ionisationEnergy = cloud_.constProps(typeIdQ).ionisationTemperature()
                                *physicoChemical::k.value();

        // calculate if an ionisation of species Q is possible
        EcQP = reactionProperties(cloud_, p, q).translationalEnergy() + reactionProperties(cloud_, p, q).EEleQ();

        if ((EcQP - ionisationEnergy) > VSMALL)
        {
            // Ionisation can occur
//             totalReactionProbability += 1.0;
//             reactionProbabilities[3] = 1.0;
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

        if (dissocReactionP)
        {
            nTotABDissociationReactions_++;
            nABDissociationReactionsPerTimeStep_++;

            if (allowSplitting_)
            {
                relax_ = false;
                
                reactionProperties(cloud_, p, q).dissociateP
                (
                    heatOfReactionDissociationAB_,
                    dissociationProducts_[0],
                    p,
                    q
                );
            }
        }

        if (dissocReactionQ)
        {
            nTotCDDissociationReactions_++;
            nCDDissociationReactionsPerTimeStep_++;

            if (allowSplitting_)
            {
                relax_ = false;
                
                reactionProperties(cloud_, p, q).dissociateQ
                (
                    heatOfReactionDissociationCD_,
                    dissociationProducts_[1],
                    p,
                    q
                );
            }
        }

        if (ionisationReactionP)
        {
            nTotABIonisationReactions_++;
            nABIonisationReactionsPerTimeStep_++;

            if (allowSplitting_)
            {
                relax_ = false;
                
                reactionProperties(cloud_, p, q).ioniseP
                (
                    heatOfReactionIonisationAB_,
                    ionisationProducts_[0],
                    p,
                    q
                );
            }
        }

        if (ionisationReactionQ)
        {
            nTotCDIonisationReactions_++;
            nCDIonisationReactionsPerTimeStep_++;

            if (allowSplitting_)
            {
                relax_ = false;
                
                reactionProperties(cloud_, p, q).ioniseQ
                (
                    heatOfReactionIonisationCD_,
                    ionisationProducts_[1],
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
        scalarList reactionProbabilities(4, 0.0);

        bool dissocReactionP = false;
        bool ionisationReactionP = false;
        bool dissocReactionQ = false;
        bool ionisationReactionQ = false;

        // 4 reactions possible
        
        // 1. Dissociation of P
        // 2. Ionisation of P
        // 3. Dissociation of Q
        // 4. Ionisations of Q

        scalar EcPQ = 0.0;
        label idP = reactionProperties(cloud_, p, q).charDissLevelP();
        label imaxP = 0;

        // calculate if a dissociation of species P is possible
        EcPQ = reactionProperties(cloud_, p, q).translationalEnergy() + reactionProperties(cloud_, p, q).EVibP();

        imaxP = EcPQ/(physicoChemical::k.value()*reactionProperties(cloud_, p, q).thetaVP());

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
        EcPQ = reactionProperties(cloud_, p, q).translationalEnergy() + reactionProperties(cloud_, p, q).EEleP();

        if ((EcPQ - ionisationEnergy) > VSMALL)
        {
//             totalReactionProbability += 1.0;
//             reactionProbabilities[1] = 1.0;
        }

        scalar EcQP = 0.0;
        label idQ = reactionProperties(cloud_, p, q).charDissLevelQ();
        label imaxQ = 0;

        // calculate if a dissociation of species Q is possible
        EcQP = reactionProperties(cloud_, p, q).translationalEnergy() + reactionProperties(cloud_, p, q).EVibQ();

        imaxQ = EcQP/(physicoChemical::k.value()*reactionProperties(cloud_, p, q).thetaVQ());

        if (imaxQ - idQ > 0)
        {
            // Dissociation can occur
            totalReactionProbability += 1.0;
            reactionProbabilities[2] = 1.0;
        }

        ionisationEnergy = cloud_.constProps(typeIdQ).ionisationTemperature()
                                    *physicoChemical::k.value();

        // calculate if an ionisation of species Q is possible
        EcQP = reactionProperties(cloud_, p, q).translationalEnergy() + reactionProperties(cloud_, p, q).EEleQ();

        if ((EcQP - ionisationEnergy) > VSMALL)
        {
//             totalReactionProbability += 1.0;
//             reactionProbabilities[3] = 1.0;
        }
        // Decide if a reaction is to occur

        if (totalReactionProbability > cloud_.rndGen().scalar01())
        {
            // A chemical reaction is to occur, choose which one

            scalarList normalisedProbabilities
                            (reactionProbabilities.size(), 0.0);
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
                            // Ionisation is to occur
                            dissocReactionQ = true;
                            break;
                        }
                        if (i == 3)
                        {
                            // Dissociation reaction is to occur
                            ionisationReactionQ = true;
                            break;
                        }
                    }
                }
            }
        }

        if (dissocReactionP)
        {
            nTotCDDissociationReactions_++;
            nCDDissociationReactionsPerTimeStep_++;

            if (allowSplitting_)
            {
                relax_ = false;
                
                reactionProperties(cloud_, p, q).dissociateP
                (
                    heatOfReactionDissociationCD_,
                    dissociationProducts_[1],
                    p,
                    q
                );
            }
        }

        if (dissocReactionQ)
        {
            nTotABDissociationReactions_++;
            nABDissociationReactionsPerTimeStep_++;

            if (allowSplitting_)
            {
                relax_ = false;
                
                reactionProperties(cloud_, p, q).dissociateQ
                (
                    heatOfReactionDissociationAB_,
                    dissociationProducts_[0],
                    p,
                    q
                );
            }
        }

        if (ionisationReactionP)
        {
            nTotCDIonisationReactions_++;
            nCDIonisationReactionsPerTimeStep_++;

            if (allowSplitting_)
            {
                relax_ = false;
                
                reactionProperties(cloud_, p, q).ioniseP
                (
                    heatOfReactionIonisationCD_,
                    ionisationProducts_[1],
                    p,
                    q
                );
            }
        }

        if (ionisationReactionQ)
        {
            nTotABIonisationReactions_++;
            nABIonisationReactionsPerTimeStep_++;

            if (allowSplitting_)
            {
                relax_ = false;
                
                reactionProperties(cloud_, p, q).ioniseQ
                (
                    heatOfReactionIonisationAB_,
                    ionisationProducts_[0],
                    p,
                    q
                );
            }
        }
    }
}


void  dissociationIonisationTypeIDissimilarSpecies::
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
        label nTotCDDissociationReactions = nTotCDDissociationReactions_;
        label nTotABIonisationReactions = nTotABIonisationReactions_;
        label nTotCDIonisationReactions = nTotCDIonisationReactions_;

        // Parallel communication
        if (Pstream::parRun())
        {
            reduce(volume, sumOp<scalar>());
            reduce(mols[0], sumOp<label>());
            reduce(mols[1], sumOp<label>());
            reduce(nTotABDissociationReactions, sumOp<label>());
            reduce(nTotCDDissociationReactions, sumOp<label>());
            reduce(nTotABIonisationReactions, sumOp<label>());
            reduce(nTotCDIonisationReactions, sumOp<label>());
        }

        numberDensities_[0] = (mols[0]*cloud().nParticle())/volume;
        numberDensities_[1] = (mols[1]*cloud().nParticle())/volume;

        const scalar deltaT = mesh_.time().deltaT().value();

        word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
        word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];

        word dissociationProductMolA =
                cloud_.typeIdList()[dissociationProducts_[0][0]];
        word dissociationProductMolB =
                cloud_.typeIdList()[dissociationProducts_[0][1]];
        word dissociationProductMolC =
                cloud_.typeIdList()[dissociationProducts_[1][0]];
        word dissociationProductMolD =
                cloud_.typeIdList()[dissociationProducts_[1][1]];

        word ionisationProductMolA =
                cloud_.typeIdList()[ionisationProducts_[0][0]];
        word ionisationProductMolB =
                cloud_.typeIdList()[ionisationProducts_[0][1]];
        word ionisationProductMolC =
                cloud_.typeIdList()[ionisationProducts_[1][0]];
        word ionisationProductMolD =
                cloud_.typeIdList()[ionisationProducts_[1][1]];

        if ((numberDensities_[0] > 0.0) && (numberDensities_[1] > 0.0))
        {
            scalar reactionRate1 = 0.0;
            scalar reactionRate2 = 0.0;
            scalar reactionRate3 = 0.0;
            scalar reactionRate4 = 0.0;

            reactionRate1 =
            (
                nTotABDissociationReactions
                * cloud_.nParticle()
            )/(counterIndex*deltaT*numberDensities_[0]
            * numberDensities_[1]*volume);

            reactionRate2 =
            (
                nTotCDDissociationReactions
                * cloud_.nParticle()
            )/(counterIndex*deltaT*numberDensities_[0]
            * numberDensities_[1]*volume);

            reactionRate3 =
            (
                nTotABIonisationReactions
                * cloud_.nParticle()
            )/(counterIndex*deltaT*numberDensities_[0]
            * numberDensities_[1]*volume);

            reactionRate4 =
            (
                nTotCDIonisationReactions
                * cloud_.nParticle()
            )/(counterIndex*deltaT*numberDensities_[0]
            * numberDensities_[1]*volume);

            Info << "Dissociation type I reaction " <<  reactantMolA
                 << " + " << reactantMolB << " --> "
                 << dissociationProductMolA << " + "
                 << dissociationProductMolB << " + " << reactantMolB
                 << ", reaction rate = " << reactionRate1
                 << nl
                 << "Dissociation type I reaction " <<  reactantMolB
                 << " + " << reactantMolA << " --> "
                 << dissociationProductMolC << " + "
                 << dissociationProductMolD << " + " << reactantMolA
                 << ", reaction rate = " << reactionRate2
                 << nl
                 << "Ionisation type I reaction " <<  reactantMolA
                 << " + " << reactantMolB << " --> "
                 << ionisationProductMolA << " + "
                 << ionisationProductMolB << " + " << reactantMolB
                 << ", reaction rate = " << reactionRate3
                 << nl
                 << "Ionisation type I reaction " <<  reactantMolB
                 << " + " << reactantMolA << " --> "
                 << ionisationProductMolC << " + "
                 << ionisationProductMolD << " + " << reactantMolA
                 << ", reaction rate = " << reactionRate4
                 << endl;
        }
    }
    else
    {
        label nTotABDissociationReactions = nTotABDissociationReactions_;
        label nTotCDDissociationReactions = nTotCDDissociationReactions_;
        label nTotABIonisationReactions = nTotABIonisationReactions_;
        label nTotCDIonisationReactions = nTotCDIonisationReactions_;

        label nABDissociationReactionsPerTimeStep =
            nABDissociationReactionsPerTimeStep_;
        label nCDDissociationReactionsPerTimeStep =
            nCDDissociationReactionsPerTimeStep_;
        label nABIonisationReactionsPerTimeStep =
            nABIonisationReactionsPerTimeStep_;
        label nCDIonisationReactionsPerTimeStep =
            nCDIonisationReactionsPerTimeStep_;

        // Parallel communication
        if (Pstream::parRun())
        {
            reduce(nTotABDissociationReactions, sumOp<label>());
            reduce(nTotCDDissociationReactions, sumOp<label>());
            reduce(nTotABIonisationReactions, sumOp<label>());
            reduce(nTotCDIonisationReactions, sumOp<label>());

            reduce(nABDissociationReactionsPerTimeStep, sumOp<label>());
            reduce(nCDDissociationReactionsPerTimeStep, sumOp<label>());
            reduce(nABIonisationReactionsPerTimeStep, sumOp<label>());
            reduce(nCDIonisationReactionsPerTimeStep, sumOp<label>());
        }

        if (nTotABDissociationReactions > VSMALL)
        {
            word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
            word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];

            word dissociationProductMolA =
                    cloud_.typeIdList()[dissociationProducts_[0][0]];
            word dissociationProductMolB =
                    cloud_.typeIdList()[dissociationProducts_[0][1]];

            Info << "Dissociation type I reaction " <<  reactantMolA << " + "
                 << reactantMolB << " --> " <<  dissociationProductMolA
                 << " + "  << dissociationProductMolB << " + "
                 << reactantMolB << " is active, nReactions this "
                 << "time step = " << nABDissociationReactionsPerTimeStep
                 << endl;
        }

        if (nTotCDDissociationReactions > VSMALL)
        {
            word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
            word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];

            word dissociationProductMolC =
                cloud_.typeIdList()[dissociationProducts_[1][0]];
            word dissociationProductMolD =
                cloud_.typeIdList()[dissociationProducts_[1][1]];

            Info << "Dissociation type I reaction " <<  reactantMolB << " + "
                 << reactantMolA << " --> " << dissociationProductMolC
                 << " + " << dissociationProductMolD << " + "
                 << reactantMolA << " is active, nReactions this "
                 << "time step = " << nCDDissociationReactionsPerTimeStep
                 << endl;
        }

        if (nTotABIonisationReactions > VSMALL)
        {
            word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
            word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];

            word ionisationProductMolA =
                cloud_.typeIdList()[ionisationProducts_[0][0]];
            word ionisationProductMolB =
                cloud_.typeIdList()[ionisationProducts_[0][1]];

            Info << "Ionisation type I reaction " <<  reactantMolA << " + "
                 << reactantMolB << " --> " << ionisationProductMolA
                 << " + " << ionisationProductMolB << " + "
                 << reactantMolB << " is active, nReactions this "
                 << "time step = " << nABIonisationReactionsPerTimeStep
                 << endl;
        }

        if (nTotCDIonisationReactions > VSMALL)
        {
            word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
            word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];

            word ionisationProductMolC =
                cloud_.typeIdList()[ionisationProducts_[1][0]];
            word ionisationProductMolD =
                cloud_.typeIdList()[ionisationProducts_[1][1]];

            Info << "Ionisation type I reaction " <<  reactantMolB << " + "
                 << reactantMolA << " --> " << ionisationProductMolC
                 << " + " << ionisationProductMolD << " + "
                 << reactantMolA << " is active, nReactions this "
                 << "time step = " << nCDIonisationReactionsPerTimeStep
                 << endl;
        }
    }

    nABDissociationReactionsPerTimeStep_ = 0;
    nCDDissociationReactionsPerTimeStep_ = 0;
    nABIonisationReactionsPerTimeStep_ = 0;
    nCDIonisationReactionsPerTimeStep_ = 0;
}


bool dissociationIonisationTypeIDissimilarSpecies::relax() const
{
    return relax_;
}

} // End namespace Foam

// ************************************************************************* //
