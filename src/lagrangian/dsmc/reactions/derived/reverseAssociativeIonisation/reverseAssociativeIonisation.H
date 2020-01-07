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
    reverseAssociativeIonisation

Description


SourceFiles
    reverseAssociativeIonisation.C


\*---------------------------------------------------------------------------*/

#ifndef reverseAssociativeIonisation_H
#define reverseAssociativeIonisation_H

#include "dsmcReaction.H"
#include "vector.H"
#include "fvMesh.H"
#include "volFields.H"
#include "dsmcCloud.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                Class reverseAssociativeIonisation Declaration
\*---------------------------------------------------------------------------*/

class reverseAssociativeIonisation
:
   public dsmcReaction
{
private:
    // Private data

        label intermediateId_;

        List<label> productIds_;
        scalar heatOfReactionDissociation_;
        scalar heatOfReactionRecombination_;
        
        label nReactions_;
        label nReactionsPerTimeStep_;

    // Private Member Functions

        void setProperties();

public:

    //- Runtime type information
    TypeName("reverseAssociativeIonisation");

    // Constructors

        //- Construct from components
        reverseAssociativeIonisation
        (
            Time& t,
            dsmcCloud& cloud,
            const dictionary& dict
        );

        // Destructor
        ~reverseAssociativeIonisation();

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
