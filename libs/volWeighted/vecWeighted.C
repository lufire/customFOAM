// Interpolation für vecWeighted bei j=-sigma*grad(vecWeighted)

#include "fvMesh.H"
#include "vecWeighted.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeSurfaceInterpolationScheme(vecWeighted)


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class Type>
void Foam::vecWeighted<Type>::clearOut()
{
    deleteDemandDrivenData(oDelta_);
    deleteDemandDrivenData(nDelta_);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const Foam::surfaceScalarField&
Foam::vecWeighted<Type>::oDelta() const	// distance face - owner
{
    if (!oDelta_)
    {
        makeDeltas();
    }

    return (*oDelta_);
}

template<class Type>
const Foam::surfaceScalarField&
Foam::vecWeighted<Type>::nDelta() const	// distance face - neighbour
{
    if (!nDelta_)
    {
        makeDeltas();
    }

    return (*nDelta_);
}

template<class Type>
void Foam::vecWeighted<Type>::makeDeltas() const // compute face-owner and face-neighbour distance
{

    const fvMesh& mesh = mesh_;

    oDelta_ = new surfaceScalarField
    (
        IOobject
        (
            "oDelta",
            mesh_.pointsInstance(),
            mesh_
        ),
        mesh_,
        dimLength
    );
    surfaceScalarField& oDelta = *oDelta_;

    nDelta_ = new surfaceScalarField
    (
        IOobject
        (
            "nDelta",
            mesh_.pointsInstance(),
            mesh_
        ),
        mesh_,
        dimLength
    );
    surfaceScalarField& nDelta = *nDelta_;

    const labelUList& owner     = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    // interpolate internal faces
    forAll(owner, facei)
    {
        oDelta[facei] = mag(mesh.C()[owner[facei]]     - mesh.Cf()[facei]);
        nDelta[facei] = mag(mesh.C()[neighbour[facei]] - mesh.Cf()[facei]);
    }

    const fvPatchList& patches = mesh.boundary();

    forAll(patches, patchi)
    {
        const fvPatch& currPatch = mesh.boundary()[patchi];
        if (currPatch.coupled()) // normal patches are just copied - delta values are needed for processor patches only
        {
            const labelUList& pOwner = mesh.boundary()[patchi].faceCells(); // owner of patch faces
            const vectorField& pCf   = mesh.Cf().boundaryField()[patchi];   // face coordinates

            forAll(pOwner, facei)
            {
                label own = pOwner[facei];

                scalar deltaGes = 1. / mesh.Cf().boundaryField()[patchi].patch().deltaCoeffs()[facei]; // Distance of cells
                oDelta.boundaryFieldRef()[patchi][facei] = mag(pCf[facei] - mesh.C()[own]); // O: owner ; a == produces an error!
                nDelta.boundaryFieldRef()[patchi][facei] = mag(deltaGes - oDelta.boundaryField()[patchi][facei]); // N: neighbour
            }
        }

        else // this part is not needed!
        {
            const labelUList& pOwner = mesh.boundary()[patchi].faceCells(); // owner of patch faces
            const vectorField& pCf   = mesh.Cf().boundaryField()[patchi];   // face coordinates

            forAll(pOwner, facei)
            {
                label own = pOwner[facei];

                oDelta.boundaryFieldRef()[patchi][facei] = mag(pCf[facei] - mesh.C()[own]); // O: owner
                nDelta.boundaryFieldRef()[patchi][facei] = 0; // N: neighbour
            }
        }
    }
}




template<class Type>
Foam::tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
Foam::vecWeighted<Type>::interpolate // interpolates cell values to faces
(
    const GeometricField<Type, fvPatchField, volMesh>& vf // volume field
) const
{
    const fvMesh& mesh = vf.mesh();
    const surfaceScalarField& oDelta = vecWeighted<Type>::oDelta();
    const surfaceScalarField& nDelta = vecWeighted<Type>::nDelta();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tvff // surface field
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
        IOobject
        (
            "vecWeighted::interpolate(" + vf.name() + ')',
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        vf.dimensions()
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& vff = tvff.ref();

    const labelUList& owner     = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    surfaceVectorField nf = mesh.Sf()/mesh.magSf();

    // interpolate internal faces
    forAll(vff, facei) // das kann evtl. beschleunigt werden; oDelta und nDelta könnte man zwischenspeichern
    {
        scalar sigmaDeltaO = (sigma_[owner[facei]] & vector(1,1,1)     & nf[facei]) / oDelta[facei];
        scalar sigmaDeltaN = (sigma_[neighbour[facei]] & vector(1,1,1) & nf[facei]) / nDelta[facei];

        vff[facei] = (vf[owner[facei]]*sigmaDeltaO + vf[neighbour[facei]]*sigmaDeltaN) / (sigmaDeltaO + sigmaDeltaN);
    }

    // interpolate patches
    typename GeometricField<Type, fvsPatchField, surfaceMesh>::Boundary& bvff = vff.boundaryFieldRef();

    forAll(bvff, patchi)
    {
        fvsPatchField<Type>& pvff = bvff[patchi];

        if (!pvff.coupled()) // if not coupled - simply copy the boundary values
        {
            pvff = vf.boundaryField()[patchi];
        }
        else // e.g. processor patches have to calculated
        {    // because they contain the neighbour field value, but the face (mean) value is needed
            const labelUList& pOwner = mesh.boundary()[patchi].faceCells(); // owner of patch faces
//        coupledFvPatchField<scalar>& sigmaPatch =
//          static_cast<coupledFvPatchField<scalar>&>(sigma_.boundaryFieldRef()[patchi]);

//        volScalarField test = fvc::average(sigma_);

            tensorField sigmaN(sigma_.boundaryField()[patchi].patchNeighbourField()); // neighbour conductivity

            Field<Type> vfO = vf.boundaryField()[patchi].patchInternalField();
            Field<Type> vfN = vf.boundaryField()[patchi].patchNeighbourField();

            forAll(pOwner, facei)
            {
            	label own = pOwner[facei];

		vector n = mesh.boundary()[patchi].Sf()[facei] / mesh.boundary()[patchi].magSf()[facei];

          	scalar sigmaDeltaO = (sigma_[own]   & vector(1,1,1) & n) / oDelta.boundaryField()[patchi][facei];
          	scalar sigmaDeltaN = (sigmaN[facei] & vector(1,1,1) & n) / nDelta.boundaryField()[patchi][facei];

          	pvff[facei] = (vfO[facei]*sigmaDeltaO + vfN[facei]*sigmaDeltaN) / (sigmaDeltaO + sigmaDeltaN);

          	//Info << nl << mesh.boundary()[patchi].name() << " " << oDelta << " "
          	//   << nDelta << " " << sigma[own] << " " << sigmaN[facei] << " " 
          	//   << vfO[facei] << " " << vfN[facei] <<  " " << deltaGes << nl << endl;
            }
        }
    }
    return tvff;
}

}
// ************************************************************************* //
