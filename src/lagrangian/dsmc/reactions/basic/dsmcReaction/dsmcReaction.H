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
    dsmcReaction

Description
    Basic/abstract class of a state controller


SourceFiles
    dsmcReactionI.H
    dsmcReaction.C
    dsmcReactionIO.C

\*---------------------------------------------------------------------------*/

#ifndef dsmcReaction_H
#define dsmcReaction_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "vector.H"
#include "volFields.H"
#include "Random.H"
#include "dsmcParcel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                        Class dsmcReaction Declaration
\*---------------------------------------------------------------------------*/

class dsmcReaction
{
protected:

    // Protected data


        const fvMesh& mesh_;

        dsmcCloud& cloud_;

        label nTotReactions_;

        label nReactionsPerTimeStep_;

        bool reactWithLists_;
        
        //- subDictionary containing the properties
        dictionary propsDict_;
        
        word reactionName_;

        List<word> reactants_;

        List<label> reactantIds_;
        
        label rDof1_;
        label rDof2_;
        label vDof1_;
        label vDof2_;
        label charge1_;
        label charge2_;
        
        bool relax_;

        // remove if not used for heat bath
        bool allowSplitting_;
        bool writeRatesToTerminal_;
        scalar volume_;
        List<scalar> numberDensities_;

      // Private Member Functions
        
        void setCommonReactionProperties();


public:
    
    //- Class to hold common properties required to perform reactions
    class reactionProperties
    {
    protected:

        // Protected data
        
        const fvMesh& mesh_;
        dsmcCloud& cloud_;


    private:

        // Private data
        
            //typeId of particle P
            label typeIdP_;
            
            //type Id of particle Q
            label typeIdQ_;

            // velocity of particle P
            vector UP_;
            
            // velocity of particle Q
            vector UQ_;
            
            // rotational energy of particle P
            scalar ERotP_;
            
             // rotational energy of particle Q
            scalar ERotQ_;
            
            // electrical energy level of particle P
            label ELevelP_;
            
             // electrical energy level of particle Q
            label ELevelQ_;

            
            // vibrational energy level of particle P
            labelList vibLevelP_;
            
             // vibrational energy level of particle Q
            labelList vibLevelQ_;
            
    public:

        // Constructors

            //- Construct from two parcels
            reactionProperties
            (
                dsmcCloud& cloud, 
                dsmcParcel& p, 
                dsmcParcel& q
            );

        // Member functions

            //- Return access to the particle P velocity
            vector& UP();
            
            //- Return access to the particle Q velocity
            vector& UQ();
            
            //- Return access to the particle P rotational energy
            scalar& ERotP();
            
            //- Return access to the particle Q rotational energy
            scalar& ERotQ();
            
            //- Return access to the particle P vibrational energy
            scalar EVibP() const;
            
            //- Return access to the particle Q vibrational energy
            scalar EVibQ() const;
            
            //- Return access to the particle P electrical energy level
            scalar EEleP() const;
            
            //- Return access to the particle Q electrical energy level
            scalar EEleQ() const;
            
            //- Return access to the particle P mass
            scalar mP() const;
            
            //- Return access to the particle Q mass
            scalar mQ() const;
            
            //- Return access to the reduced mass
            scalar mR() const;
            
            //- Return access to the translational energy
            scalar translationalEnergy() const;
            
            //- Return access to omegaPQ
            scalar omegaPQ() const;
            
            //- Return access to particle P thetaD
            scalar thetaDP() const;
            
            //- Return access to particle Q thetaD
            scalar thetaDQ() const;
            
            //- Return access to particle P thetaV
            scalar thetaVP() const;
            
            //- Return access to particle Q thetaV
            scalar thetaVQ() const;
            
            //- Return access to particle P Zref
            scalar ZrefP() const;
            
            //- Return access to particle Q Zref
            scalar ZrefQ() const;
            
            //- Return access to particle P refTempZv
            scalar refTempZvP() const;
            
            //- Return access to particle Q refTempZv
            scalar refTempZvQ() const;
            
            //- Return access to species P dissociation level
            label charDissLevelP() const;
            
            //- Return access to species Q dissociation level
            label charDissLevelQ() const;
            
            //- Return access to species P number of electronic energy levels
            label jMaxP() const;
            
            //- Return access to species Q number of electronic energy levels
            label jMaxQ() const;
            
            //- Return access to species P rotational energy modes
            label rotationalDofP() const;
            
            //- Return access to species Q rotational energy modes
            label rotationalDofQ() const;
            
            // Return centre of mass velocity of particle pair
            vector Ucm() const;
            
            //- Return access to species P electronic energies
            const List<scalar>& EEListP() const;
            
            //- Return access to species Q electronic energies
            const List<scalar>& EEListQ() const;
            
            //- Return access to species P degeneracies
            const List<label>& gListP() const;
            
            //- Return access to species Q degeneracies
            const List<label>& gListQ() const;
            
            //- Dissociate particle P
            void dissociateP
            (
                const scalar& heatOfReaction,
                const List<label>& productIds,
                dsmcParcel& p, 
                dsmcParcel& q
            );
            
            //- Dissociate particle Q
            void dissociateQ
            (
                const scalar& heatOfReaction,
                const List<label>& productIds,
                dsmcParcel& p, 
                dsmcParcel& q
            );
            
            //- Ionise particle P
            void ioniseP
            (
                const scalar& heatOfReaction,
                const List<label>& productIds,
                dsmcParcel& p, 
                dsmcParcel& q
            );
            
            //- Ionise particle Q
            void ioniseQ
            (
                const scalar& heatOfReaction,
                const List<label>& productIds,
                dsmcParcel& p, 
                dsmcParcel& q
            );
            
            //- Perform associative ionisation
            void associativeIonisation
            (
                const scalar& heatOfReactionIntermediateIonisation,
                const scalar& heatOfReactionRecombination,
                const List<label>& assIonProductIds,
                dsmcParcel& p, 
                dsmcParcel& q
            );
    };

    //- Runtime type information
    TypeName("dsmcReaction");

    // Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            dsmcReaction,
            dictionary,
            (
                Time& t,
                dsmcCloud& cloud,
                const dictionary& dict
            ),
            (t, cloud, dict)
        );

    // Constructors

        //- Construct from components
        dsmcReaction
        (
            Time& t,
            dsmcCloud& cloud,
            const dictionary& dict
        );


    // Selectors

        static autoPtr<dsmcReaction> New
        (
            Time& t,
            dsmcCloud& cloud,
            const dictionary& dict
        );


    // Destructor

        virtual ~dsmcReaction();

    // Member Functions

        virtual void initialConfiguration() = 0;

        virtual void reaction
        (
            dsmcParcel& p,
            dsmcParcel& q
        ) = 0;

        virtual void reaction
        (
            dsmcParcel& p,
            dsmcParcel& q,
            const DynamicList<label>& candidateList,
            const List<DynamicList<label> >& candidateSubList,
            const label candidateP,
            const List<label>& whichSubCell
        ) = 0;

        virtual bool relax() const = 0;

        virtual bool tryReactMolecules(const label typeIdP,
                                       const label typeIdQ) = 0;

        virtual void outputResults(const label counterIndex) = 0;

        // Access
        const dsmcCloud& cloud() const;

        label nTotReactions() const;

        label nReactionsPerTimeStep() const;
        label& nReactionsPerTimeStep();

        bool reactWithLists() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "dsmcReactionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
