    fileName path(UMean.rootPath()/UMean.caseName()/"graphs"/UMean.instance());
    mkDir(path);

    vectorField UMeanValues(channelIndexing.collapse(UMean));
    scalarField productionValues(channelIndexing.collapse(production));
    scalarField dissipationValues(channelIndexing.collapse(dissipation));
    scalarField vpgValues(channelIndexing.collapse(vpg));
    scalarField molDiffValues(channelIndexing.collapse(molDiff));
    scalarField turbDiffValues(channelIndexing.collapse(turbDiff));
    

    const scalarField& y = channelIndexing.y();

    makeGraph(y, UMeanValues.component(0), "Ux", path, gFormat); 
    makeGraph(y, productionValues, "P_ii", path, gFormat); 
    makeGraph(y, dissipationValues, "E_ii", path, gFormat); 
    makeGraph(y, vpgValues, "Pi_ii", path, gFormat); 
    makeGraph(y, molDiffValues, "D_ii", path, gFormat); 
    makeGraph(y, turbDiffValues, "T_ii", path, gFormat); 
