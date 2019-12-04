// Interpolation für volWeighted bei j=-sigma*grad(volWeighted)

#include "fvMesh.H"
#include "volWeighted.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

makeSurfaceInterpolationScheme(volWeighted)

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class Type>
void Foam::volWeighted<Type>::clearOut()
{
    deleteDemandDrivenData(oDelta_);
    deleteDemandDrivenData(nDelta_);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const Foam::surfaceScalarField&
Foam::volWeighted<Type>::oDelta() const	// distance face - owner
{
    if (!oDelta_)
    {
        makeDeltas();
    }

    return (*oDelta_);
}

template<class Type>
const Foam::surfaceScalarField&
Foam::volWeighted<Type>::nDelta() const	// distance face - neighbour
{
    if (!nDelta_)
    {
        makeDeltas();
    }

    return (*nDelta_);
}

template<class Type>
void Foam::volWeighted<Type>::makeDeltas() const // compute face-owner and face-neighbour distance
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

    const surfaceVectorField n = mesh.Sf()/mesh.magSf();

    // interpolate internal faces
    // all distances are NORMAL to the face
    forAll(owner, facei)
    {
        oDelta[facei] = mag(n[facei] & (mesh.C()[owner[facei]]     - mesh.Cf()[facei]));
        nDelta[facei] = mag(n[facei] & (mesh.C()[neighbour[facei]] - mesh.Cf()[facei]));
    }

    const fvPatchList& patches = mesh.boundary();

    forAll(patches, patchi)
    {
        const fvPatch& currPatch = mesh.boundary()[patchi];
        const vectorField np = currPatch.Sf()/currPatch.magSf(); // patch normal vector

        if (currPatch.coupled()) // normal patches are just copied - delta values are needed for processor patches only
        {
            const labelUList& pOwner = mesh.boundary()[patchi].faceCells(); // owner of patch faces
            const vectorField& pCf   = mesh.Cf().boundaryField()[patchi];   // face coordinates

            forAll(pOwner, facei)
            {
                label own = pOwner[facei];

		// all distances are NORMAL to the face!
                oDelta.boundaryFieldRef()[patchi][facei] = mag(np[facei] & (pCf[facei] - mesh.C()[own])); // O: owner ; a == produces an error!
            }
            // weight: w = delta_neighbour / delta in ORTHOGONAL direction
            nDelta.boundaryFieldRef()[patchi] = currPatch.weights() * oDelta.boundaryFieldRef()[patchi] / (1.-currPatch.weights());
        }
        else // this part is not needed!
        {
            const labelUList& pOwner = mesh.boundary()[patchi].faceCells(); // owner of patch faces
            const vectorField& pCf   = mesh.Cf().boundaryField()[patchi];   // face coordinates

            forAll(pOwner, facei)
            {
                label own = pOwner[facei];
		// all distances are NORMAL to the face!
                oDelta.boundaryFieldRef()[patchi][facei] = mag(np[facei] & (pCf[facei] - mesh.C()[own])); // O: owner
                nDelta.boundaryFieldRef()[patchi][facei] = mag(np[facei] & (pCf[facei] - mesh.C()[own])); // N: neighbour
            }
        }
    }
}




template<class Type>
Foam::tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
Foam::volWeighted<Type>::interpolate // interpolates cell values to faces
(
    const GeometricField<Type, fvPatchField, volMesh>& vf // volume field
) const
{
    const fvMesh& mesh = vf.mesh();
    const surfaceScalarField& oDelta = volWeighted<Type>::oDelta();
    const surfaceScalarField& nDelta = volWeighted<Type>::nDelta();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tvff // surface field
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "volWeighted::interpolate(" + vf.name() + ')',
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

    // interpolate internal faces
    forAll(vff, facei) 
    {
        scalar sigmaDeltaO = sigma_[owner[facei]]     / oDelta[facei];
        scalar sigmaDeltaN = sigma_[neighbour[facei]] / nDelta[facei];

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
            scalarField sigmaN(sigma_.boundaryField()[patchi].patchNeighbourField()); // neighbour conductivity

            Field<Type> vfO = vf.boundaryField()[patchi].patchInternalField();
            Field<Type> vfN = vf.boundaryField()[patchi].patchNeighbourField();

            forAll(pOwner, facei)
            {
                label own = pOwner[facei];

                scalar sigmaDeltaO = sigma_[own]   / oDelta.boundaryField()[patchi][facei];
                scalar sigmaDeltaN = sigmaN[facei] / nDelta.boundaryField()[patchi][facei];

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
