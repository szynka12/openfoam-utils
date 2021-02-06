    fileName path(UMean.rootPath()/UMean.caseName()/"graphs"/UMean.instance());
    mkDir(path);

    scalarField y = channelIndexing.y();
    vectorField UMeanValues(channelIndexing.collapse(UMean));
    
    // compute U_tau
    scalar u_tau = Foam::sqrt(nu * UMeanValues[0].component(0) / y[0] );

    // normalise U
    UMeanValues = UMeanValues / u_tau;
    y = y * u_tau / nu;
    
    scalarField productionValues(channelIndexing.collapse(production));
    scalarField dissipationValues(channelIndexing.collapse(dissipation));
    scalarField vpgValues(channelIndexing.collapse(vpg));
    scalarField molDiffValues(channelIndexing.collapse(molDiff));
    scalarField turbDiffValues(channelIndexing.collapse(turbDiff));
    


    makeGraph(y, UMeanValues.component(0), "Ux", path, gFormat); 
    makeGraph(y, productionValues, "P_ii", path, gFormat); 
    makeGraph(y, dissipationValues, "E_ii", path, gFormat); 
    makeGraph(y, vpgValues, "Pi_ii", path, gFormat); 
    makeGraph(y, molDiffValues, "D_ii", path, gFormat); 
    makeGraph(y, turbDiffValues, "T_ii", path, gFormat); 
