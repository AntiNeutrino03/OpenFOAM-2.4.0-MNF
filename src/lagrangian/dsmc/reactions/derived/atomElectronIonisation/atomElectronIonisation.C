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

#include "atomElectronIonisation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(atomElectronIonisation, 0);

addToRunTimeSelectionTable
(
    dsmcReaction,
    atomElectronIonisation,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

atomElectronIonisation::atomElectronIonisation
(
    Time& t,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcReaction(t, cloud, dict),
    productIdsIon_(),
    heatOfReactionIon_(),
    nIonisationReactions_(0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

atomElectronIonisation::~atomElectronIonisation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void atomElectronIonisation::initialConfiguration()
{
    setCommonReactionProperties();
    setProperties();
}


void atomElectronIonisation::setProperties()
{
    if (reactants_[0] == reactants_[1])
    {
        FatalErrorIn("atomElectronIonisation::setProperties()")
            << "Reactant molecules cannot be same species." << nl
            << exit(FatalError);
    }

    // check that reactant one is an 'ATOM'

    const label rDof1 =
            cloud_.constProps(reactantIds_[0]).rotationalDegreesOfFreedom();
            
    const label vDof1 =
        cloud_.constProps(reactantIds_[0]).vibrationalDegreesOfFreedom();

    if (rDof1 > VSMALL)
    {
        FatalErrorIn("atomElectronIonisation::setProperties()")
            << "First reactant must be an atom "
            << "(not a molecule or an electron): " << reactants_[0]
            << nl
            << exit(FatalError);
    }

    if (vDof1 > VSMALL)
    {
         FatalErrorIn("atomElectronIonisation::setProperties()")
            << "Reactions are currently only implemented "
            << "for monatomic and diatomic species."
            << " This is a polyatomic:" << reactants_[0]
            << nl
            << exit(FatalError);
    }
    
    const label charge1 = cloud_.constProps(reactantIds_[0]).charge();

    if (charge1 == -1)
    {
        FatalErrorIn("atomElectronIonisation::setProperties()")
            << "First reactant must not be an electron "
            << reactants_[1]
            << nl
            << exit(FatalError);
    }

    // check that reactant two is an 'ELECTRON'

    const label charge2 = cloud_.constProps(reactantIds_[1]).charge();

    if (charge2 != -1)
    {
        FatalErrorIn("atomElectronIonisation::setProperties()")
            << "Second reactant must be an electron, not "
            << reactants_[1]
            << nl
            << exit(FatalError);
    }

    // reading in ionisation products

    List<word> productMoleculesIonisation
        (propsDict_.lookup("productsOfIonisedAtom"));

    if (productMoleculesIonisation.size() != 2)
    {
        FatalErrorIn("atomElectronIonisation::setProperties()")
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
            FatalErrorIn("atomElectronIonisation::setProperties()")
                << "There should be two products (for the ionising molecule "
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
                FatalErrorIn("atomElectronIonisation::setProperties()")
                    << "Cannot find type id: "
                    << productMoleculesIonisation[r]
                    << nl
                    << exit(FatalError);
            }
        }

        // check that product one is a 'ATOM'

        const scalar rDof2 =
            cloud_.constProps(productIdsIon_[0]).rotationalDegreesOfFreedom();
            
        const scalar vDof2 =
            cloud_.constProps(productIdsIon_[0]).vibrationalDegreesOfFreedom();

        if (rDof2 > 1)
        {
            FatalErrorIn("atomElectronIonisation::setProperties()")
                << "First product must be an atom (not an atom): "
                << productMoleculesIonisation[0]
                << nl
                << exit(FatalError);
        }
        
        if (vDof2 > VSMALL)
        {
            FatalErrorIn("atomElectronIonisation::setProperties()")
                << "Reactions are currently only implemented "
                << "for monatomic and diatomic species."
                << " This is a polyatomic:" << productMoleculesIonisation[0]
                << nl
                << exit(FatalError);
        }

        // check that product two is an 'ELECTRON'

        const label charge = cloud_.constProps(productIdsIon_[1]).charge();

        if (charge != -1)
        {
            FatalErrorIn("atomElectronIonisation::setProperties()")
                << "Second product must be an electron: "
                << productMoleculesIonisation[1]
                << nl
                << exit(FatalError);
        }
    }
    
    heatOfReactionIon_ =
                readScalar(propsDict_.lookup("heatOfReactionIonisation"));
}


bool atomElectronIonisation::tryReactMolecules
(
    label typeIdP,
    const label typeIdQ
)
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


void atomElectronIonisation::reaction
(
    dsmcParcel& p,
    dsmcParcel& q,
    const DynamicList<label>& candidateList,
    const List<DynamicList<label> >& candidateSubList,
    const label candidateP,
    const List<label>& whichSubCell
)
{}


void atomElectronIonisation::reaction
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

        scalar ionisationEnergy =
            cloud_.constProps(typeIdP).ionisationTemperature()
           *physicoChemical::k.value();

        scalar EcPPIon = reactionProperties(cloud_, p, q).translationalEnergy() + reactionProperties(cloud_, p, q).EEleP();

        if ((EcPPIon - ionisationEnergy) > VSMALL)
        {
            nIonisationReactions_++;
            nReactionsPerTimeStep_++;

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

        scalar ionisationEnergy =
            cloud_.constProps(typeIdQ).ionisationTemperature()
            *physicoChemical::k.value();

        scalar EcPQIon = reactionProperties(cloud_, p, q).translationalEnergy() + reactionProperties(cloud_, p, q).EEleQ();

        if ((EcPQIon - ionisationEnergy) > VSMALL)
        {
             nIonisationReactions_++;
             nReactionsPerTimeStep_++;

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


void  atomElectronIonisation::outputResults(const label counterIndex)
{
    if (writeRatesToTerminal_ == true)
    {
        const List<DynamicList<dsmcParcel*> >& cellOccupancy =
            cloud_.cellOccupancy();

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
        label nTotReactionsIonisation = nIonisationReactions_;

        // Parallel communication
        if (Pstream::parRun())
        {
            reduce(volume, sumOp<scalar>());
            reduce(mols[0], sumOp<label>());
            reduce(mols[1], sumOp<label>());
            reduce(nTotReactionsIonisation, sumOp<label>());
        }

        numberDensities_[0] = (mols[0]*cloud().nParticle())/volume;
        numberDensities_[1] = (mols[1]*cloud().nParticle())/volume;

        const scalar deltaT = mesh_.time().deltaT().value();

        const word& reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
        const word& reactantMolB = cloud_.typeIdList()[reactantIds_[1]];

        const word& productMolC = cloud_.typeIdList()[productIdsIon_[0]];
        const word& productMolD = cloud_.typeIdList()[productIdsIon_[1]];

        if ((numberDensities_[0] > 0.0) && (numberDensities_[1] > 0.0))
        {
            scalar reactionRateIonisation = 0.0;

            reactionRateIonisation =
            (
                nTotReactionsIonisation
               *cloud_.nParticle()
            )
           /(
                counterIndex*deltaT*numberDensities_[0]
               *numberDensities_[1]*volume
            );

            Info<< "Electron ionisation reaction "
                <<  reactantMolA << " + " << reactantMolB
                <<  " --> "
                << productMolC << " + " << productMolD << " + "
                << reactantMolB
                << ", reaction rate = " << reactionRateIonisation
                << endl;
        }
    }
    else
    {
        label nTotReactionsIonisation = nIonisationReactions_;
        label nReactionsPerTimeStep = nReactionsPerTimeStep_;

        if (Pstream::parRun())
        {
            reduce(nTotReactionsIonisation, sumOp<label>());
            reduce(nReactionsPerTimeStep, sumOp<label>());
        }

        if (nTotReactionsIonisation > VSMALL)
        {
            const word& reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
            const word& reactantMolB = cloud_.typeIdList()[reactantIds_[1]];

            const word& productMolA = cloud_.typeIdList()[productIdsIon_[0]];
            const word& productMolB = cloud_.typeIdList()[productIdsIon_[1]];

            Info<< "Electron ionisation reaction "
                <<  reactantMolA << " + " << reactantMolB
                <<  " --> "
                << productMolA << " + " << productMolB << " + "
                << reactantMolB
                << " is active, nReactions this time step = "
                << nReactionsPerTimeStep << endl;
        }
    }

    nReactionsPerTimeStep_ = 0.0;
}


bool atomElectronIonisation::relax() const
{
    return relax_;
}

} // End namespace Foam

// ************************************************************************* //
