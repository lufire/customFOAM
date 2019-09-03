{
    Info<< "Solving chemistry" << endl;

    chemistry.solve
    (
        runTime.deltaT().value()
    );


    // turbulent time scale
    if (turbulentReaction)
    {
        volScalarField tk =
                Cmix*sqrt(turbulence->nuEff()/turbulence->epsilon());
        volScalarField tc = chemistry.tc();

        // Chalmers PaSR model
        gamma = (runTime.deltaT() + tc)/(runTime.deltaT() + tc + tk);
    }
    else
    {
        gamma = 1.0;
    }

    chemistrySh = gamma*chemistry.Sh()();
}
