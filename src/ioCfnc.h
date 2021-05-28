#ifndef ioCfnc__h
#define ioCfnc__h
// This file is for io_fnc that require ioClasses, too
#include "ioClass.h"
#include "io_fnc.h"

std::tuple<ioPads, TH2D*, TH1D*, TH1D*, TH1D*> ioMakeClosure(
        const char* file, 
        const char* rooResponseName,
        ioOptMap opt={},
        ioOptMap dict={{
            /* "which_truth","A",   // other options are B and all */
            /* "which_response","B", */
            "iter_unfold",3,
            "StatsResponse",0,
            "TitleResponse","",
            "name_A","", // unfold A with B for closure
            "name_B","",
            "TitleResponse","",
            "xTitleResponse","#it{p}_{T} Measured",
            "yTitleResponse","#it{p}_{T} Truth",
            "xTitleRatio","",
            "yTitleTruth","#frac{1}{N_{jets}} #frac{dN_{jets}}{d#it{p}_{T}}",
            "TitleTruth","black:truth red:unfold",
            "yTitleRatio","Ratio: unf./truth"}}
        );


#endif
