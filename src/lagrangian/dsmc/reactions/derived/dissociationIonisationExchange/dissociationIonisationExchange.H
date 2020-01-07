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

Class
    dissociationIonisationExchange

Description


SourceFiles
    dissociationIonisationExchange.C


\*---------------------------------------------------------------------------*/

#ifndef dissociationIonisationExchange_H
#define dissociationIonisationExchange_H

#include "dsmcReaction.H"
#include "vector.H"
#include "fvMesh.H"
#include "volFields.H"
#include "dsmcCloud.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
               Class dissociationIonisationExchange Declaration
\*---------------------------------------------------------------------------*/

class dissociationIonisationExchange
:
   public dsmcReaction
{
private:
    // Private data

        List<label> exchangeProductIds_;
        List<label> chargeExchangeProductIds_;
        List<label> dissociationProductIds_;
        List<label> ionisationProductsIdsP_;
        List<label> ionisationProductsIdsQ_;

        bool chargedAtom_;
        bool chargedMolecule_;
        bool chargeExchange_;

        scalar heatOfReactionDiss_;
        scalar heatOfReactionExch_;
        scalar heatOfReactionChargeExchange_;
        scalar heatOfReactionIonP_;
        scalar heatOfReactionIonQ_;

        scalar aCoeff_;
        scalar bCoeff_;
        scalar aCoeffCharge_;
        scalar bCoeffCharge_;

        label nTotExchangeReactions_;
        label nTotChargeExchangeReactions_;
        label nTotDissociationReactions_;
        label nTotIonisationReactionsP_;
        label nTotIonisationReactionsQ_;

        label nExchangeReactionsPerTimeStep_;
        label nChargeExchangeReactionsPerTimeStep_;
        label nDissociationReactionsPerTimeStep_;
        label nIonisationReactionsPPerTimeStep_;
        label nIonisationReactionsQPerTimeStep_;

    // Private Member Functions

        void reactExchangeMolecule
        (
            dsmcParcel& p,
            label newTypeId,
            const label newEVibLevel,
            const scalar newERot,
            const vector& newU
        );

        void reactExchangeAtom
        (
            dsmcParcel& p,
            label newTypeId,
            const vector& newU
        );

        void setProperties();


public:

    //- Runtime type information
    TypeName("dissociationIonisationExchange");

    // Constructors

        //- Construct from components
        dissociationIonisationExchange
        (
            Time& t,
            dsmcCloud& cloud,
            const dictionary& dict
        );

        // Destructor
        ~dissociationIonisationExchange();

    // Member Functions

        // initial configuration - read in properties from chemReactDict
        void initialConfiguration();

        // apply a reaction between parcel p and q
        void reaction
        (
            dsmcParcel& p,
            dsmcParcel& q,
            const DynamicList<label>& candidateList,
            const List<DynamicList<label> >& candidateSubList,
            const label candidateP,
            const List<label>& whichSubCell
        );

        void reaction
        (
            dsmcParcel& p,
            dsmcParcel& q
        );

        // return bool if molecules belong to this reaction model
        bool tryReactMolecules(label typeIdP, const label typeIdQ);

        // Access

        // return bool for relax = true (run normal collision model)
        // relax = false (apply reaction to molecules)
        bool relax() const;

        void outputResults(const label counterIndex);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
