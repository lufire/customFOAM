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
//const volScalarField& psi = thermo.psi();

volScalarField& hs = thermo.he();  

volScalarField& T = thermo.T();

#include "createPhi.H"
//#include "compressibleCreatePhi.H"

//===============================
// particle interaction modelling
//===============================

Info<< "\nReading momentum exchange field Ksl\n" << endl;
volScalarField Ksl
(
    IOobject
    (
        "Ksl",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
    //dimensionedScalar("0", dimensionSet(1, -3, -1, 0, 0), 1.0)
);

Info<< "\nReading voidfraction field voidfraction = (Vgas/Vparticle)\n" << endl;
volScalarField voidfraction
(
    IOobject
    (
        "voidfraction",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "\nReading solid velocity field field Us\n" << endl;
volVectorField Us
(
    IOobject
    (
        "Us",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField m
(
    IOobject
    (
        "m",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh
);

//===============================

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
label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);

singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);
// Set the turbulence into the reaction model
//reaction->setTurbulence(turbulence());

Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt
(
    IOobject
    (
        "dpdt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("dpdt", p.dimensions()/dimTime, 0)
);

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));

multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(Y, i)
{
    fields.add(Y[i]);
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

//// Relative permittivity [F/m]
//dimensionedScalar epsilonR
//(
//    transportProperties.lookup("epsilonR")
//);
//
//// Vacuum permittivity [C]
//const dimensionedScalar epsilon0
//(
//    "epsilon0",
//    dimensionSet(-1, -3, 4, 0, 0, 2, 0),
//    8.854187817e-12
//);  

const dimensionedScalar& R = constant::physicoChemical::R;
const dimensionedScalar& F = constant::physicoChemical::F;
const dimensionedScalar& kB = constant::physicoChemical::k;
const dimensionedScalar& e = constant::electromagnetic::e;

// create cfdemCloud
myCfdemCloud particleCloud(mesh);
        
// Create electrochemistry model
electrochemistryModels::electrochemistryModel electrochemistry(thermo, particleCloud);
PtrList<electrochemistryModels::electrodeModels::electrodeModel>& 
    electrodes(electrochemistry.electrodes());
electrochemistryModels::electrolyteModels::electrolyteModel& 
    electrolyte(electrochemistry.electrolyte());

//volScalarField iTrans
//(
//    IOobject
//    (
//        "iTrans",
//        runTime.timeName(),
//        mesh,
//        IOobject::NO_READ,
//        IOobject::AUTO_WRITE
//    ),
//    mesh,
//    dimensionedScalar("iTrans", dimCurrent/dimVolume, 0.0)
//);
//volScalarField& phiEl = electrochemistry.phiEl();

volScalarField& phiEl = electrolyte.phiE();
volVectorField& id = electrolyte.iD();
volVectorField& il = electrolyte.i();
//volScalarField phiEl
//(
//    IOobject
//    (
//        "phiEl",
//        runTime.timeName(),
//        mesh,
//        IOobject::MUST_READ,
//        IOobject::AUTO_WRITE
//    ),
//    mesh
//);
//
//volScalarField phiEs
//(
//    IOobject
//    (
//        "phiEs",
//        runTime.timeName(),
//        mesh,
//        IOobject::MUST_READ,
//        IOobject::AUTO_WRITE
//    ),
//    mesh
//);


//volVectorField il
//(
//    IOobject
//    (
//        "il",
//        runTime.timeName(),
//        mesh,
//        IOobject::NO_READ,
//        IOobject::AUTO_WRITE
//    ),
//    mesh,
//    dimensionedVector
//    (
//        "il",
//        dimensionSet(0,-2,0,0,0,1,0),
//        vector(0.0, 0.0, 0.0)
//    )
//);
//
//volVectorField id
//(
//    IOobject
//    (
//        "id",
//        runTime.timeName(),
//        mesh,
//        IOobject::NO_READ,
//        IOobject::AUTO_WRITE
//    ),
//    mesh,
//    dimensionedVector
//    (
//        "id",
//        dimensionSet(0,-2,0,0,0,1,0),
//        vector(0.0, 0.0, 0.0)
//    )
//);

//volVectorField is
//(
//    IOobject
//    (
//        "is",
//        runTime.timeName(),
//        mesh,
//        IOobject::NO_READ,
//        IOobject::AUTO_WRITE
//    ),
//    mesh,
//    dimensionedVector
//    (
//        "is",
//        dimensionSet(0,-2,0,0,0,1,0),
//        vector(0.0, 0.0, 0.0)
//    )
//);

//volVectorField im
//(
//    IOobject
//    (
//        "im",
//        runTime.timeName(),
//        mesh,
//        IOobject::NO_READ,
//        IOobject::AUTO_WRITE
//    ),
//    mesh,
//    dimensionedVector
//    (
//        "im",
//        dimensionSet(0,-2,0,0,0,1,0),
//        vector(0.0, 0.0, 0.0)
//    )
//);

//volScalarField sigma
//(
//    IOobject
//    (
//        "sigma",
//        runTime.timeName(),
//        mesh,
//        IOobject::NO_READ,
//        IOobject::AUTO_WRITE
//    ),
//    mesh,
//    dimensionedScalar
//    (
//        "sigma", 
//        sqr(dimCurrent)*pow3(dimTime)/(dimVol*dimMass),
//        1000.0
//    )
//);

//volScalarField kappa
//(
//    IOobject
//    (
//        "kappa",
//        runTime.timeName(),
//        mesh,
//        IOobject::NO_READ,
//        IOobject::AUTO_WRITE
//    ),
//    mesh,
//    dimensionedScalar
//    (
//        "kappa", 
//        sqr(dimCurrent)*pow3(dimTime)/(dimVol*dimMass),
//        50.0
//    )
//);


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

//volScalarField div_il
//(
//    IOobject
//    (
//        "div_il",
//        runTime.timeName(),
//        mesh,
//        IOobject::NO_READ,
//        IOobject::AUTO_WRITE
//    ),
//    mesh,
//    dimensionedScalar
//    (
//        "div_il",
//        dimensionSet(0,-3,0,0,0,1,0),
//        0.0 
//    )
//);
//
//volScalarField div_is
//(
//    IOobject
//    (
//        "div_is",
//        runTime.timeName(),
//        mesh,
//        IOobject::NO_READ,
//        IOobject::AUTO_WRITE
//    ),
//    mesh,
//    dimensionedScalar
//    (
//        "div_is",
//        dimensionSet(0,-3,0,0,0,1,0),
//        0.0 
//    )
//);
// ------------------------------------------------------------------------- //
//autoPtr<liquidMultiSpeciesTransportModel> mstm
//(
//    liquidMultiSpeciesTransportModel::New
//    (
//        thermo,
//        turbulence()
//    )
//);
//mstm().updateCoefficients();

// Initialize species charge number [-], Molar mass [kg/mol],
// concentration [mol/m³], diffusion coefficients [m²/s], 
// and mobility [m²mol/J-s]

//Info << "C size " << Y.size() << endl;
//const PtrList<label>& z = mstm().z();
//PtrList<dimensionedScalar> M(C.size());
//PtrList<volVectorField> N;
//N.setSize(Y.size());
//const PtrList<volScalarField>& D = mstm().D();
//PtrList<volScalarField> DEff;
//DEff.setSize(D.size());
////Info << "D size " << D.size() << endl;
//PtrList<volScalarField> mu(Y.size()-1);

const volScalarField& kappa(electrolyte.kappa());
const PtrList<volScalarField>& D(electrolyte.D());
const PtrList<label>& z = electrolyte.z();
const PtrList<volScalarField>& t = electrolyte.transNumberModel().t();
PtrList<dimensionedScalar> M(C.size());
PtrList<volVectorField>& N = electrolyte.N();
PtrList<volScalarField> mu(D.size());
//PtrList<volScalarField> Deff(D.size());
//PtrList<volScalarField> tEff(t.size());

scalar BruggemanCoeff = 1.5;
volScalarField porosityFactor = pow(voidfraction, BruggemanCoeff);

volScalarField kappaEff = kappa*porosityFactor;
// Species molar weight in kg/mol
forAll(Y, i)
{
    M.set
    (
        i,
        new dimensionedScalar
        (
           "M_"+Y[i].name(),
           dimensionSet(dimMass/dimMoles),
           thermo.composition().W(i)/1000.0
        )
    );
    mu.set
    (
        i,
        new volScalarField(D[i]/(R*T))
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
    F*z[0]*D[0]
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
    D[0]
);
////label nChargedSpecies = 0;
//forAll(Y, i)
//{
//    if (z[i] != 0)
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
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zC", dimMoles/dimVol, 0.0)
);

scalar iAvgTime = 0.0;
scalar phiEsAvgTime = 0.0;
label timeStepCounter = 0;
label startCounter = 50;

const IOdictionary& electrochemDict = electrochemistry.dict();
const word transportSolver = electrochemDict.lookup("transportSolver");
