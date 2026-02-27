#!/usr/bin/env python

# Parametrise a TH1D or TH1F histo with Chebychev polynomials with
# a fit of the coefficients

from ROOT import TFile, TCanvas, TGraphErrors, gStyle, TF1, gInterpreter
from ROOT import TH1D, TH1F, TLegend
from ROOT import gROOT

gROOT.SetBatch(True) # Enable batch mode (no graphics windows)

# Chebychev from C++ class
gInterpreter.ProcessLine( '#include "ChebychevApproximation.hh"' )
from ROOT import ChebychevApproximation

# Get histos from file
def getHistosFromFile( keys, files ):
    hists= dict()
    for f in files:
        tf= TFile.Open( f )
        for key in keys:
            hist= tf.Get( key )
            if hist:
                hist.SetDirectory( 0 )
                if len(files) > 1:
                    hists[f.split(".")[0]+"_"+key]= hist
                else:
                    hists[key]= hist
            else:
                raise RuntimeError( "getHistosFromFile: histogram for "+key+" not found" )
        tf.Close()
    return hists

# Study parametrisations
def mtopDependence( mlbhists, n=9, a=40.0, b=160.0, txt="Parametrisation2", opt="nf" ):

    lnorm= False
    if "n" in opt:
        lnorm= True        

    mlbhistKeys= sorted( mlbhists.keys() )
    key= mlbhistKeys[0]
    hist= mlbhists[key]
    binw= hist.GetBinWidth( 1 )
    nbins= int((b-a)/binw)
    print( "mtopDependence: bins in range", a, "to", b, "from", key, ":", nbins )
    
    # Do the parametrisation of the hists
    global fittfs, chebtfs, cas
    fittfs= dict()
    chebtfs= dict()
    cas= list()
    fitParameters= dict()
    fitParErrors= dict()
    for key in mlbhistKeys:
        mlbhist= mlbhists[key]
        print( "mtopDependence: histo ", key )
        ca= ChebychevApproximation( mlbhist, a, b, n, lnorm )
        cas.append( ca )
        pars= ca.getCoefficients()
        fittf= TF1( key+"_tf1", ca, a, b, n )
        chebtf= TF1( key+"_cheb_tf1", ca, a, b, n )
        fittf.Print()
        fittfs[key]= fittf        
        chebtfs[key]= chebtf
        for i in range( n ):            
            fittf.SetParameter( i, pars[i] )
            fittf.SetParName( i, "c"+str(i) )
            chebtf.SetParameter( i, pars[i] )
            chebtf.SetParName( i, "c"+str(i) )
        fitresult= mlbhist.Fit( fittf, "0S", "", a, b )
        cov= fitresult.GetCorrelationMatrix()
        cov.Print()
        pars= fittf.GetParameters()
        errs= fittf.GetParErrors()            
        fitpars= list()
        fiterrs= list()
        for i in range( n ):
            fitpars.append( pars[i] )
            fiterrs.append( errs[i] )
        fitParameters[key]= fitpars
        fitParErrors[key]= fiterrs

    # Results of parametrisation
    print( "mtopDependence: Results for fitted Chebychev coefficients" )
    for key in mlbhistKeys:
        pars= fitParameters[key]
        errs= fitParErrors[key]
        print( key, end= ": " )
        for i in range( n ):
            value= pars[i]
            error= errs[i]
            print( "{0:7.3f} +/- {1:7.3f}".format( value, error ), end=" " )
        print()
    if not lnorm:
        print( "Calculated normalised parameters:" )
        for key in mlbhistKeys:
            pars= fitParameters[key]
            errs= fitParErrors[key]
            norm= pars[0]
            for i in range( n ):
                value= pars[i]/norm
                error= errs[i]/norm
                print( "{0:7.4f} {1:7.4f}".format( value, error ), end=" " )
            print()
        
    # Some globals so plots stay alive
    global canv1, canv2, tges, ratioPlots

    # Fit info on plots
    gStyle.SetOptFit( 1111 )
    from ROOT import kBlue
    
    # Plots of parametrisations
    canv1= TCanvas( "canv1", "Fit plots", 600, 800 )    
    canv1.Divide( 2, 3 )
    icanv= 0
    for key in mlbhistKeys:
        icanv+= 1
        canv1.cd( icanv )
        mlbhist= mlbhists[key]
        mlbhist.SetMarkerStyle(20)
        mlbhist.SetMarkerSize(0.5)
        mlbhist.Draw("P")
        chebtf= chebtfs[key]
        chebtf.SetLineColor( kBlue )
        chebtf.Draw( "same" )
        fittf= fittfs[key]
        fittf.Draw( "same" )
        leg1= TLegend( 0.3, 0.7, 0.5, 0.85 )
        leg1.AddEntry(mlbhist, "Data", "P")
        leg1.AddEntry(chebtf, "Chebyshev (pre)", "L")
        leg1.AddEntry(fittf, "Chebyshev (post)", "L")
        leg1.SetBorderSize(0)
        leg1.DrawClone("same")
        #leg1.Clear()
    del leg1
    pdffilename= txt
    if lnorm:
        pdffilename+= "_N"
    pdffilename+= "_"+str(n)+"_"+str(a)+"_"+str(b)+".pdf"
    canv1.Print( pdffilename+"(" )
    
    # Ratio plots
    canv1.SetTitle( "Significance plots" )
    icanv= 0
    ratioPlots= list()
    for key in mlbhistKeys:
        icanv+= 1
        pad= canv1.cd( icanv )
        pad.SetGrid()
        mlbhist= mlbhists[key]
        fittf= fittfs[key]
        ax= mlbhist.GetXaxis()
        firstbin= ax.FindBin( a )
        lastbin= ax.FindBin( b )
        ratioPlot= TGraphErrors()
        ratioPlot.SetName( key+"_ratio" )
        ratioPlot.SetTitle( "Significance "+key )
        ratioPlots.append( ratioPlot )
        for ibin in range( firstbin, lastbin+1 ):
            value= mlbhist.GetBinContent( ibin )
            error= mlbhist.GetBinError( ibin )
            bincenter= mlbhist.GetBinCenter( ibin )
            fitvalue= fittf.Eval( bincenter )
            significance= (fitvalue-value)/error
            ipoint= ibin-firstbin
            ratioPlot.SetPoint( ipoint, bincenter, significance )
            ratioPlot.SetPointError( ipoint, 0.0, 1.0 )
        ratioPlot.SetMinimum( -3.0 )
        ratioPlot.SetMaximum( 4.0 )
        ratioPlot.GetXaxis().SetLimits( a, b )
        ratioPlot.SetMarkerStyle( 20 )
        ratioPlot.SetMarkerSize( 0.75 )
        ratioPlot.Fit( "pol0" )
        ratioPlot.Draw( "ap" )
    canv1.Print( pdffilename )
    
    # Plots of coefficients vs mtop
    tges= list()
    for icoeff in range( n ):
        tge= TGraphErrors()
        tge.SetName( "tge"+str(icoeff) )
        tge.SetTitle( "Chebychev coeff c"+str(icoeff) )
        imtopPoint= 0
        for key in mlbhistKeys:
            if "WbWb" in key:
                mtop= float( key.split("_")[3] ) # Get top mass from file name
            else:
                tokens= key.split("_")
                mtop= float( tokens[len(tokens)-1] )
            fitpars= fitParameters[key]
            fiterrs= fitParErrors[key]
            value= fitpars[icoeff]
            error= fiterrs[icoeff]
            tge.SetPoint( imtopPoint, mtop, value )
            tge.SetPointError( imtopPoint, 0.0, error )
            imtopPoint+= 1
        tges.append( tge )
    canv2= TCanvas( "canv2", "Coefficients mtop dependence", 800, 800 )
    canv2.SetTitle( "Coefficients mtop dependence" )
    canv2.DivideSquare( n )
    for icoeff in range( n ):
        canv2.cd( icoeff+1 )
        tge= tges[icoeff]
        tge.SetMarkerColor( 1 )
        tge.SetMarkerStyle( 20 )
        tge.SetMarkerSize( 0.75 )
        tge.Fit( "pol1" )
        tge.Draw( "ap" )
    canv2.Print( pdffilename+")" )
                
    return

def parametrisationPlots( opt="" ):
    
    mlbhistKeys= [ "h_mlb_171.0", "h_mlb_172.0", "h_mlb_172.5", "h_mlb_173.0", "h_mlb_174.0" ]
    mlbhists= getHistosFromFile( mlbhistKeys, ["output_signal.root"] )
    mtopDependence( mlbhists, n=8, a=40.0, b=148.0, txt="ParametrisationRoot_ttbar", opt=opt )
    mtopDependence( mlbhists, n=8, a=40.0, b=160.0, txt="ParametrisationRoot_ttbar", opt=opt )
    mtopDependence( mlbhists, n=9, a=40.0, b=148.0, txt="ParametrisationRoot_ttbar", opt=opt )
    mtopDependence( mlbhists, n=9, a=40.0, b=160.0, txt="ParametrisationRoot_ttbar", opt=opt )
    mtopDependence( mlbhists, n=10, a=40.0, b=148.0, txt="ParametrisationRoot_ttbar", opt=opt )
    mtopDependence( mlbhists, n=10, a=40.0, b=160.0, txt="ParametrisationRoot_ttbar", opt=opt )

    # Histos Broken at the moment
    # mlbhistKeys= [ "tReco_mlb_minavg_171.0", "tReco_mlb_minavg_172.0",
    #                "tReco_mlb_minavg_172.5", "tReco_mlb_minavg_173.0", 
    #                "tReco_mlb_minavg_174.0" ]
    # mlbhists= getHistosFromFile( mlbhistKeys, "output_signal.root" )
    # mtopDependence( mlbhists, n=8, a=40.0, b=148.0, txt="Parametrisation2_WWbb", opt=opt )
    # mtopDependence( mlbhists, n=8, a=40.0, b=160.0, txt="Parametrisation2_WWbb", opt=opt )
    # mtopDependence( mlbhists, n=9, a=40.0, b=148.0, txt="Parametrisation2_WWbb", opt=opt )
    # mtopDependence( mlbhists, n=9, a=40.0, b=160.0, txt="Parametrisation2_WWbb", opt=opt )
    
    return

def parametrisationPlotsWbWb( opt="" ):

    mlbhistKeys= [ "m_bl_fine" ]
    files = ["WbWb_Slurm_Template_165.0.root", "WbWb_Slurm_Template_167.5.root", "WbWb_Slurm_Template_170.0.root", "WbWb_Slurm_Template_172.5.root", "WbWb_Slurm_Template_175.0.root"]
    mlbhists= getHistosFromFile( mlbhistKeys, files )
    mtopDependence( mlbhists, n=8, a=50.0, b=180.0, txt="ParametrisationRoot_WbWb", opt=opt )
    mtopDependence( mlbhists, n=8, a=50.0, b=300.0, txt="ParametrisationRoot_WbWb", opt=opt )
    mtopDependence( mlbhists, n=9, a=50.0, b=180.0, txt="ParametrisationRoot_WbWb", opt=opt )
    mtopDependence( mlbhists, n=9, a=50.0, b=300.0, txt="ParametrisationRoot_WbWb", opt=opt )
    mtopDependence( mlbhists, n=10, a=50.0, b=180.0, txt="ParametrisationRoot_WbWb", opt=opt )
    mtopDependence( mlbhists, n=10, a=50.0, b=300.0, txt="ParametrisationRoot_WbWb", opt=opt )

def main():
    parametrisationPlotsWbWb( opt="n" )
    return


# Run as main
if __name__ == '__main__':
   main()

   
