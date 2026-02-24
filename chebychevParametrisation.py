#!/usr/bin/env python

# Parametrise a TH1D or TH1F histo with Chebychev polynomials with
# a fit of the coefficients

from ROOT import TFile, TCanvas, TGraphErrors, gStyle, TF1, gInterpreter
from ROOT import TH1D, TH1F

# Chebychev from C++ class
gInterpreter.ProcessLine( '#include "ChebychevApproximation.hh"' )
from ROOT import ChebychevApproximation

# Get histos from file
def getHistosFromFile( keys, filename ):
    tf= TFile.Open( filename )
    hists= dict()
    for key in keys:
        hist= tf.Get( key )
        if hist:
            hist.SetDirectory( 0 )
            hists[key]= hist
        else:
            raise RuntimeError( "getHistosFromFile: histogram for "+key+" not found" )
    tf.Close()
    return hists

# Study parametrisations
def mtopDependence( mlbhists, n=9, a=40.0, b=160.0, txt="Parametrisation2", opt="n" ):

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
    global fittfs, cas
    fittfs= dict()
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
        fittf.Print()
        fittfs[key]= fittf        
        for i in range( n ):            
            fittf.SetParameter( i, pars[i] )
            fittf.SetParName( i, "c"+str(i) )
        mlbhist.Fit( fittf, "0", "", a, b )
        pars= fittf.GetParameters()
        errs= fittf.GetParErrors()
        fitpars= list()
        fiterrs= list()
        for i in range( n ):
            fitpars.append( pars[i] )
            fiterrs.append( errs[i] )
        fitParameters[key]= fitpars
        fitParErrors[key]= fiterrs
        
    print( "mtopDependence: Results for Chebychev coefficients" )
    for key in mlbhistKeys:
        pars= fitParameters[key]
        errs= fitParErrors[key]
        for i in range( n ):
            value= pars[i]
            error= errs[i]
            print( "{0:7.4f} {1:7.4f}".format( value, error ), end=" " )
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

    # Plots of parametrisations
    canv1= TCanvas( "canv1", "Fit plots", 600, 800 )    
    canv1.Divide( 2, 3 )
    icanv= 0
    for key in mlbhistKeys:
        icanv+= 1
        canv1.cd( icanv )
        mlbhist= mlbhists[key]
        mlbhist.Draw()
        fittf= fittfs[key]
        fittf.Draw( "same" )

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
    mlbhists= getHistosFromFile( mlbhistKeys, "output_signal.root" )
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


def main():
    parametrisationPlots( opt="n" )
    return


# Run as main
if __name__ == '__main__':
   main()

   
