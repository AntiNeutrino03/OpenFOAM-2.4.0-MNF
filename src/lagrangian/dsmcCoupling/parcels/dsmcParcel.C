/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
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

#include "dsmcParcel.H"
#include "dsmcCloud.H"
#include "meshTools.H"
#include "dsmcVolFields.C"

//#include "solidParticleCloud.H"
//#include "solidParticle.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTemplateTypeNameAndDebug(Cloud<solidParticle>, 0);
    defineTemplateTypeNameAndDebug(Cloud<dsmcParcel>, 0); 
    //defination of the solidParticleCouplingCloud class, which is a templeted cloud of solid particles.
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dsmcParcel::move
(
    dsmcParcel::trackingData& td,
    const scalar trackTime
)
{
    td.switchProcessor = false;
    td.keepParticle = true;

    if(stuckToWall() == 0)
    {
        const polyMesh& mesh = td.cloud().pMesh();
        const polyBoundaryMesh& pbMesh = mesh.boundaryMesh();

        if(newParcel() == 1)
        {
            Random& rndGen(td.cloud().rndGen());
            stepFraction() = rndGen.scalar01(); 
            newParcel() = 0;
        }
        
        
        scalar tEnd = (1.0 - stepFraction())*trackTime;
        const scalar dtMax = tEnd;
                    
        // For reduced-D cases, the velocity used to track needs to be
        // constrained, but the actual U_ of the parcel must not be
        // altered or used, as it is altered by patch interactions an
        // needs to retain its 3D value for collision purposes.
        vector Utracking = U_;
            
        while (td.keepParticle && !td.switchProcessor && tEnd > ROOTVSMALL)
        {           
            Utracking = U_;

            if(!td.cloud().axisymmetric())
            {
                // Apply correction to position for reduced-D cases, 
                // but not axisymmetric cases
                meshTools::constrainToMeshCentre(mesh, position());
                
                // Apply correction to velocity to constrain tracking for
                // reduced-D cases,  but not axisymmetric cases
                meshTools::constrainDirection
                (
                    mesh, mesh.solutionD(), Utracking
                );
            }

            // Set the Lagrangian time-step
            scalar dt = min(dtMax, tEnd);
            
            //dt *= rayTrace(position() + dt*Utracking, td);
            dt *= trackToFace(position() + dt*Utracking, td);
            
            tEnd -= dt;

            stepFraction() = 1.0 - tEnd/trackTime;
                
            // - face tracking info
            if( face() != -1 )    //*******
            {
                //--  measure flux properties
                td.cloud().tracker().updateFields
                (
                    *this
                );
            }

            if (onBoundary() && td.keepParticle)
            {
                if (isA<processorPolyPatch>(pbMesh[patch(face())]))
                {
                    td.switchProcessor = true;
                }

                forAll(td.cloud().boundaries().cyclicBoundaryModels(), c)
                {
                    const labelList& faces =td.cloud().boundaries().             
                                    cyclicBoundaryModels()[c]->allFaces();

                    if(findIndex(faces, this->face()) != -1)
                    {
                        td.cloud().boundaries().
                            cyclicBoundaryModels()[c]->controlMol(*this, td);
                    }
                }
            }
        }
    }

    return td.keepParticle;
}

bool Foam::dsmcParcel::hitPatch
(
    const polyPatch&,
    trackingData& td,
    const label,
    const scalar,
    const tetIndices&
)
{
    return false;
}

bool Foam::dsmcParcel::hitPatch
(
    const polyPatch&,
    trackingData& td,
    const label
)
{
    return false;
}


void Foam::dsmcParcel::hitProcessorPatch
(
    const processorPolyPatch&,
    trackingData& td
)
{
    td.switchProcessor = true;
}

void Foam::dsmcParcel::hitWallPatch
(
    const wallPolyPatch& wpp,
    trackingData& td,
    const tetIndices& tetIs
)
{
    //-find which patch has been hit
    label patchIndex = wpp.index();

    const label& patchModelId = td.cloud().boundaries().
    patchToModelIds()[patchIndex];

    // apply a boundary model when a particle collides with this poly patch
    td.cloud().boundaries().
    patchBoundaryModels()[patchModelId]->controlParticle(*this, td);
}

void Foam::dsmcParcel::hitWallPatch
(
    const wallPolyPatch& wpp,
    trackingData& td
)
{
//     Info << "PERFORMING WALL CODE" << endl;
    //-find which patch has been hit
    label patchIndex = wpp.index();

    const label& patchModelId = td.cloud().boundaries().
    patchToModelIds()[patchIndex];

    // apply a boundary model when a particle collides with this poly patch
    td.cloud().boundaries().
    patchBoundaryModels()[patchModelId]->controlParticle(*this, td);
}

void Foam::dsmcParcel::hitPatch
(
    const polyPatch& pp,
    trackingData& td
)
{
    //-find which patch has been hit
    label patchIndex = pp.index();

    const label& patchModelId = td.cloud().boundaries().
    patchToModelIds()[patchIndex];

    // apply a boundary model when a particle collides with this poly patch
    td.cloud().boundaries().
    patchBoundaryModels()[patchModelId]->controlParticle(*this, td);
}


 //template<class ParcelType>
void Foam::dsmcParcel::transformProperties
(
    const tensor& T
)
{
   particle::transformProperties(T);
   U_ = transform(T, U_);
}


 //template<class ParcelType>
void Foam::dsmcParcel::transformProperties
(
    const vector& separation
)
{
  particle::transformProperties(separation);
}

bool Foam::solidParticle::move 
//move function includes the implementation of the Reynolds number, needed for the calculation of the drag force and the new velocity. 
(
    trackingData& td,
    const scalar trackTime
)
{
    td.switchProcessor = false;
    td.keepParticle = true;

    const polyBoundaryMesh& pbMesh = mesh_.boundaryMesh();

    scalar tEnd = (1.0 - stepFraction())*trackTime;
    scalar dtMax = tEnd;

    while (td.keepParticle && !td.switchProcessor && tEnd > SMALL)
    {
        if (debug)
        {
            Info<< "Time = " << mesh_.time().timeName()
                << " trackTime = " << trackTime
                << " tEnd = " << tEnd
                << " steptFraction() = " << stepFraction() << endl;
        }

        // set the lagrangian time-step
        scalar dt = min(dtMax, tEnd);

        // remember which cell the parcel is in
        // since this will change if a face is hit
        //label cellI = cell();

        dt *= trackToFace(position() + dt*Up_, td);

        tEnd -= dt;
        stepFraction() = 1.0 - tEnd/trackTime;

        //cellPointWeight cpw(mesh_, position(), cellI, face());
        //scalar rhoc = td.rhoInterp().interpolate(cpw);
        //rhoc --> the carrier phase density.
        //vector Uc = td.UInterp().interpolate(cpw);
        //Uc --> the carrier phase velocity.
        //scalar nuc = td.nuInterp().interpolate(cpw);
        //nuc --> the carrier phase viscosity.

//        scalar rhop = td.cloud().rhop();


//         bool Foam::dsmcParcel::move
//         (
//             dsmcParcel::trackingData& td,
//             const scalar trackTime
//         )
        
//         vector Uc = U_; //dsmc parcel velocity
// 
//         scalar magUr = mag(Uc - Up_); //relative velocity
/*
        scalar ReFunc = 1.0;
        scalar Re = magUr*d_/nuc;

        if (Re > 0.01)
        {
            ReFunc += 0.15*pow(Re, 0.687);
        }

        scalar Dc = (24.0*nuc/d_)*ReFunc*(3.0/4.0)*(rhoc/(d_*rhop));

        U_ = (U_ + dt*(Dc*Uc + (1.0 - rhoc/rhop)*td.g()))/(1.0 + dt*Dc); //velocity update
*/




// { 
//     sampleCounter_++;
// 
//     const scalar& nParticle = cloud_.nParticle();
// //     rhoNInstantaneous_ = 0.0;
// 
//     if(sampleInterval_ <= sampleCounter_)
//     {
//         nTimeSteps_ += 1.0;
// 
// //         if(densityOnly_)
// //         {
//             forAllConstIter(dsmcCloud, cloud_, iter)
//             {
//                 const dsmcParcel& p = iter();
//                 label iD = findIndex(typeIds_, p.typeId());
// 
//                 if(iD != -1)
//                 {
//                     const label& cell = p.cell();
//                     const scalar& cellVolume = mesh_.cellVolumes()[cell];
//                     const scalar& nAvTimeSteps = nTimeSteps_;
//                     const scalar& mass = cloud_.constProps(p.typeId()).mass();
// 
//                     rhoNMean_[cell] += 1.0;
// //                     rhoNInstantaneous_[cell] += 1.0;
// 
// //                     if(cloud_.axisymmetric())
// //                     {
// //                         const point& cC = cloud_.mesh().cellCentres()[cell];
// // 
// //                         //scalar radius = cC.y();
// //                         scalar radius = sqrt((p.position().y()*p.position().y()) + (p.position().z()*p.position().z()));
// //                         scalar RWF = 1.0;
// //                         RWF = 1.0 + 
// //                             cloud_.maxRWF()*(radius/cloud_.radialExtent());
// // 
// //                         rhoNMeanXnParticle_[cell] += (RWF*nParticle);
// //                         rhoMMeanXnParticle_[cell] += (mass*RWF*nParticle);
// //                     }
// //                     else
// //                     {
//                         rhoNMeanXnParticle_[cell] += nParticle;
//                         rhoMMeanXnParticle_[cell] += (mass*nParticle);
// //                     }
//                         rhoN_[cell] = 
//                         (rhoNMeanXnParticle_[cell])/(nAvTimeSteps*cellVolume);
//                 }
//             }
// //         }
//     }
// }


//         const scalar& mass = cloud_.constProps(p.typeId()).mass();
//         
//         forAll(rhoN_, c)
//         {
//                         scalar F_ = (mass*rhoN_[c])*(3.14*pow((d_/2),2))*magUr*(((1+(4/9)*(1-epsilon_)*(1-alpha_))*abs(magUr))+((1-epsilon_)*alpha_*((pow(3.14,0.5))/3)*(pow((2*1.380658e-23*Tp_)/mass_),2)));
//                 
//         }



        

        if (onBoundary() && td.keepParticle)
        {
            if (isA<processorPolyPatch>(pbMesh[patch(face())]))
            {
                td.switchProcessor = true;
            }
        }
    }

    return td.keepParticle;
}


bool Foam::solidParticle::hitPatch 
//some additional functions are implemented that determine what happens when a patch is hit by a particle and these are different depending on whether it is a processor patch or a wall patch for example.
(
    const polyPatch&,
    trackingData&,
    const label,
    const scalar,
    const tetIndices&
)
{
    return false;
}


void Foam::solidParticle::hitProcessorPatch
(
    const processorPolyPatch&,
    trackingData& td
)
{
    td.switchProcessor = true;
}


void Foam::solidParticle::hitWallPatch
(
    const wallPolyPatch& wpp,
    trackingData& td,
    const tetIndices& tetIs
)
{
    vector nw = tetIs.faceTri(mesh_).normal();
    nw /= mag(nw);

    scalar Un = Up_ & nw;
    vector Ut = Up_ - Un*nw;

    if (Un > 0)
    {
//        Up_ -= (1.0 + td.cloud().e())*Un*nw;
    }

//    Up_ -= td.cloud().mu()*Ut;
}


void Foam::solidParticle::hitPatch
(
    const polyPatch&,
    trackingData& td
)
{
    td.keepParticle = false;
}


void Foam::solidParticle::transformProperties (const tensor& T)
{
    particle::transformProperties(T);
    Up_ = transform(T, Up_);
}


void Foam::solidParticle::transformProperties(const vector& separation)
{
    particle::transformProperties(separation);
}


Foam::scalar Foam::solidParticle::wallImpactDistance(const vector&) const
{
//    return 0.5*d_;
}

// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

#include "dsmcParcelIO.C"
//#include "solidParticleIO.C"

// ************************************************************************* //
