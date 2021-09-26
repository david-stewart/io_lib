#ifndef ioCfnc__h
#define ioCfnc__h

// This file is for io_fnc that require ioClasses, too
#include "ioOptMap.h"
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

RooUnfoldResponse ioMakeRooUnfoldResponse(
        const char* name, const char* file, const char* tag_M, const char* tag_T );

TH1D* ioBlankTH1D(ioBinVec bins={{0.,1.}}, ioOptMap options={}, bool draw_vert_lines=false);


#endif
