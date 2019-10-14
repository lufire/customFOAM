//Info<< nl << "Creating reaction model" << endl;
//
//autoPtr<combustionModels::concCombustionModel> reaction // OF 2.3
//(
//    combustionModels::concCombustionModel::New(mesh) // OF 2.3
//);

autoPtr<concChemistryModel> pChemistry // OF 2.3
(
    concChemistryModel::New(mesh) // OF 2.3
);

concChemistryModel& chemistry = pChemistry(); // OF 2.3

concReactionThermo& thermo = chemistry.thermo(); // OF 2.3
thermo.validate(args.executable(), "h", "e");

concBasicMultiComponentMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();
PtrList<volScalarField>& C = composition.C();

word inertSpecie(thermo.lookup("inertSpecie"));

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()//*dimensionedScalar("1",dimMass/dimVol,1.0)
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField& p = thermo.p();

volScalarField& hs = thermo.he();  

volScalarField& T = thermo.T();

#include "createPhi.H"
//#include "compressibleCreatePhi.H"

//Info << "Creating turbulence model.\n" << nl;
//autoPtr<compressible::turbulenceModel> turbulence
//(
//    compressible::turbulenceModel::New
//    (
//        rho,
//        U,
//        phi,
//        thermo
//    )
//);

// Set the turbulence into the reaction model
//reaction->setTurbulence(turbulence());

singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

//Info<< "Creating field dpdt\n" << endl;
//volScalarField dpdt
//(
//    IOobject
//    (
//        "dpdt",
//        runTime.timeName(),
//        mesh
//    ),
//    mesh,
//    dimensionedScalar("dpdt", p.dimensions()/dimTime, 0)
//);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("PIMPLE"), pRefCell, pRefValue);

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));

multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(Y, speciesI)
{
    fields.add(Y[speciesI]);
}
fields.add(hs);

DimensionedField<scalar, volMesh> chemistrySh
(
    IOobject
    (
        "chemistry::Sh",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("chemistrySh", dimEnergy/dimTime/dimVolume, 0.0)
);

volScalarField dQ
(
    IOobject
    (
        "dQ",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("dQ", dimEnergy/dimTime, 0.0)
);

volScalarField gamma 
(
    IOobject
    (
        "gamma",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);

//IOdictionary transportProperties
//(
//    IOobject
//    (
//        "transportProperties",
//        runTime.constant(),
//        mesh,
//        IOobject::MUST_READ,
//        IOobject::NO_WRITE
//    )
//);



const dimensionedScalar& R = constant::physicoChemical::R;
const dimensionedScalar& F = constant::physicoChemical::F;
const dimensionedScalar& kB = constant::physicoChemical::k;
const dimensionedScalar& e = constant::electromagnetic::e;


// Create electrochemistry model
electrochemistryModels::electrochemistryModel electrochemistry(thermo);
PtrList<electrochemistryModels::electrodeModels::electrodeModel>& 
    electrodes(electrochemistry.electrodes());

electrochemistryModels::electrolyteModels::electrolyteModel& 
    electrolyte(electrochemistry.electrolyte());

volScalarField div_il 
(
    IOobject
    (

        "div_il",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("div_il", dimCurrent/dimVolume, 0.0)
);
//volScalarField iTrans
//(
//    IOobject
//    (
//
//        "iTrans",
//        runTime.timeName(),
//        mesh,
//        IOobject::NO_READ,
//        IOobject::AUTO_WRITE
//    ),
//    mesh,
//    dimensionedScalar("iTrans", dimCurrent/dimVolume, 0.0)
//);
volScalarField& phiEl = electrolyte.phiE();
//volVectorField gradPhiE = fvc::grad(phiEl);

volVectorField& id = electrolyte.iD();
volVectorField& il = electrolyte.i();


volVectorField E
(
    IOobject
    (
        "E",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "E",
        dimensionSet(1,1,-3,0,0,-1,0),
        vector(0.0, 0.0, 0.0)
    )
);

// Initialize species charge number [-], Molar mass [kg/mol],
// concentration [mol/m³], diffusion coefficients [m²/s], 
// and mobility [m²mol/J-s]

const volScalarField& kappa(electrolyte.kappa());
const PtrList<volScalarField>& D(electrolyte.D());
const PtrList<label>& z = electrolyte.z();
const PtrList<volScalarField>& t = electrolyte.transNumberModel().t();
PtrList<dimensionedScalar> M(C.size());
PtrList<volVectorField>& N = electrolyte.N();
PtrList<volScalarField> Deff(D.size());
PtrList<volScalarField> mu(D.size());

// Species molar weight in kg/mol
forAll(Y, speciesI)
{
    M.set
    (
        speciesI,
        new dimensionedScalar
        (
           "M_"+Y[speciesI].name(),
           dimensionSet(dimMass/dimMoles),
           thermo.composition().W(speciesI)/1000.0
        )
    );
    Deff.set
    (
        speciesI,
        new volScalarField 
        (
            IOobject
            (
                "Deff_"+Y[speciesI].name(),
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            D[speciesI]
        )
    );
    mu.set
    (
        speciesI,
        new volScalarField(Deff[speciesI]/(R*T))
    );
}

volScalarField Dm 
(
    IOobject
    (
        "Dm",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    F*z[0]*Deff[0]
);

volScalarField Di
(
    IOobject
    (
        "Di",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    Deff[0]
);

//label nChargedSpecies = 0;
//forAll(Y, speciesI)
//{
//    if (z[speciesI] != 0)
//    {
//        nChargedSpecies++;
//    }
//}

volScalarField rhoE
(
    IOobject
    (
        "rhoE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("rhoE", dimCurrent*dimTime/dimVol, 0.0)
);

volScalarField zC 
(
    IOobject
    (
        "zC",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zC", dimMoles/dimVol, 0.0)
);

//// Relative permittivity [F/m]
//const dimensionedScalar& epsilonR = electrolyte.epsilon();
//
//// Vacuum permittivity [C]
//const dimensionedScalar& epsilon0 = constant::electromagnetic::epsilon0;

scalar iAvgTime = 0.0;
scalar phiEsAvgTime = 0.0;
label timeStepCounter = 0;
label startCounter = 50;

const IOdictionary& electrochemDict = electrochemistry.dict();
const word transportSolver = electrochemDict.lookup("transportSolver");

//// Dimensionless quantities
//
//dimensionedScalar L_ref("L_ref", dimLength, 1e-3);
//dimensionedScalar T_ref("T_ref", dimTemperature, T[0]);
//dimensionedScalar u_ref("u_ref", dimLength/dimTime, 1.0);
//dimensionedScalar t_ref = L_ref/u_ref;
//dimensionedScalar phiE_ref = kB/e*T_ref;
//const dimensionedScalar& C_ref = electrolyte.Cstd();
////dimensionedScalar kappa_ref = F*C_ref/tau_ref*sqr(L_ref)/phiE_ref;
//
//PtrList<volScalarField> C_star(C.size());
//forAll(C, speciesI)
//{
//    C_star.set
//    (
//        speciesI,
//        new volScalarField 
//        (
//            IOobject
//            (
//                "C_star_"+Y[speciesI].name(),
//                runTime.timeName(),
//                mesh,
//                IOobject::NO_READ,
//                IOobject::NO_WRITE
//            ),
//            mesh,
//            dimensionedScalar("C_star", dimless, 0.0)
//        )
//    );
//    C_star[speciesI] = C[speciesI]/C_ref;
//}
//PtrList<volVectorField> N_star(N.size());
//forAll(N, speciesI)
//{
//    N_star.set
//    (
//        speciesI,
//        new volVectorField 
//        (
//            IOobject
//            (
//                "N_star_"+Y[speciesI].name(),
//                runTime.timeName(),
//                mesh,
//                IOobject::NO_READ,
//                IOobject::NO_WRITE
//            ),
//            mesh,
//            dimensionedScalar("N_star", dimLength/dimTime, 0.0)
//        )
//    );
//}
//
//volScalarField phiE_star = phiEl/phiE_ref;
//
//volVectorField i_star = il/C_ref;
//volVectorField id_star = id/C_ref;
 
