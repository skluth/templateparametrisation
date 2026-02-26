#!/usr/bin/env python

# Parametrise a TH1 histo with Chebychev polynomials, optionally with
# a fit of coefficients

from ROOT import RooRealVar, RooGaussian, RooChebychev, RooArgList, RooFit, RooAddPdf, RooArgSet, RooHistPdf, RooAbsReal, RooLinkedList, RooCmdArg, RooDataHist, RooGenericPdf

from ROOT import kBlue, kRed, kGreen

from ROOT.RooFit import Label, Format, AutoPrecision, LineColor, Name

from ROOT import TFile, TCanvas, TGraph, TGraphErrors, gStyle, TH1D, TMath, TRatioPlot, TLegend

from ROOT import SetOwnership

from math import fabs, sin


# Test Chebyshev classes with wrapper function (object)
from ROOT import TF1
def chebyshevApprox( n=8, fopt="R" ):
    
    # gInterpreter.ProcessLine( '#include "SinFuncObj.hh"' )
    # from ROOT import SinFuncObj
    # Try sin function
    # sfo= SinFuncObj()
    def sfo( x ):
        return sin( x )
    ca= ChebychevApproximation( sfo, 1.0, 2.0, 3 )
    print( sin(1.5), ca.evaluate(1.5) )
    
    # Function around TH1D histo
    mlbhistKeys= [ "h_mlb_172.5" ]
    mlbhists= getHistosFromFile( mlbhistKeys, "output_signal.root" )
    mlbhist= mlbhists[ "h_mlb_172.5" ]
    def mlbfunc( x ):
        return mlbhist.Interpolate( x )
    mlbca= ChebychevApproximation( mlbfunc, 40.0, 160.0, n )
    print( mlbfunc( 100.0 ), mlbca.evaluate( 100.0 ) )
    pars= mlbca.getCoefficients()
    print( mlbca( 100.0, pars ) )
    # Setup TF1 and fit histo
    mlbcatf= TF1( "mlbcatf", mlbca, 40.0, 160.0, n )
    for i in range( n ):
        mlbcatf.SetParameter( i, pars[i] )
        parname= "c"+str(i)
        mlbcatf.SetParName( i, parname )
    mlbcatf.Print()
    print( mlbcatf.Eval( 100.0 ) )
    mlbhist.Fit( mlbcatf, fopt )
        
    return mlbca

# First 8 polynomials for tests of recursive polynomials
def T0( x ):
    return 1.0
def T1( x ):
    return x
def T2( x ):
    return 2.0*x**2-1.0
def T3( x ):
    return 4.0*x**3-3.0*x
def T4( x ):
    return 8.0*x**4-8.0*x**2+1.0
def T5( x ):
    return 16.0*x**5-20.0*x**3+5.0*x
def T6( x ):
    return 32.0*x**6-48.0*x**4+18.0*x**2-1.0
def T7( x ):
    return 64*x**7-112.0*x**5+56*x**3-7.0*x

# Chebychev approximation a la http://mathworld.wolfram.com/ChebyshevApproximationFormula.html
# Recursive polynomials function
def chbchvpoly( x, n ):
    result= None
    if n == 0:
        result= 1.0
    elif n == 1:
        result= x
    else:
        tnminusone= 1.0
        tn= x
        for i in range( 2, n+1 ):
            tnplusone= 2.0*x*tn - tnminusone
            tnminusone= tn
            tn= tnplusone
        result= tnplusone
    return result
# Transformations [-1,1] <-> [a,b]
def tTox( t, a, b ):
    return (b-a)/2.0*t + (a+b)/2.0
def xTot( x, a, b ):
    return (x-(a+b)/2.0)/(b-a)*2.0
# Class initialises Chebychev coefficients for given function, range and order
class ChebychevApproximation:

    def __init__( self, func, a, b, n ):
        self.__a= a
        self.__b= b
        self.__n= n
        from math import cos
        # Interpolation nodes on [-1,1] and [a,b]
        tnodes= []
        xnodes= []
        for k in range( 1, n+1 ):
            tnode= cos( (k-0.5)/n*3.14159 )
            xnode= tTox( tnode, a, b )
            tnodes.append( tnode )
            xnodes.append( xnode )
        # Coefficients c0 to cn:
        coeffs= []
        for j in range( n ):
            coeff= 0.0
            for k in range( n ):
                xnode= xnodes[k]
                tnode= tnodes[k]
                coeff+= func( xnode )*chbchvpoly( tnode, j )
            coeff*= 2.0/n
            coeffs.append( coeff )
        coeffs[0]/= 2.0
        self.__coeffs= coeffs
        print( "ChebychevApproximation: degree "+str(n)+", coefficients" )
        print( coeffs )
        # Check on nodes
        print( "ChebychevApproximation: degree "+str(n)+", approximation on nodes" )
        for k in range( n ):
            xnode= xnodes[k]
            print( xnode, func( xnode ), self.evaluate( xnode ) )
        return

    def getCoefficients( self ):
        return list( self.__coeffs )

    def evaluate( self, x ):
        return self.__call__( [ x ], self.__coeffs )
    
    # Callable with external parameters for TF1, we get plain cppyy buffers
    # as call arguments
    def __call__( self, x, pars ):
        xx= x[0]
        value= 0.0
        for j in range( self.__n ):
            value+= pars[j]*chbchvpoly( xTot( xx, self.__a, self.__b ), j )
        return value


# Tests with ChebychevApproximation
def chbchv_test( n=2, a=-5.0, b=5.0 ):
    
    from math import cos

    # Test function for approximation
    def fun( x ):
        return TMath.Gaus( x, 0.0, 2.0 )
        # return cos( x )

    # Nodes on [-1,1] and [a,b]
    tnodes= []
    xnodes= []
    for k in range( 1, n+1 ):
        tnode= cos( (k-0.5)/n*3.14159 )
        xnode= (b-a)/2.0*tnode + (a+b)/2.0
        tnodes.append( tnode )
        xnodes.append( xnode )
    print( "tnodes", tnodes )
    print( "xnodes", xnodes )
    print( "function", [ fun(xnode) for xnode in xnodes ] )

    # c0 to cn:
    coeffs= []
    for j in range( 0, n ):
        coeff= 0.0
        for k in range( n ):
            xnode= xnodes[k]
            tnode= tnodes[k]
            coeff+= fun( xnode )*chbchvpoly( tnode, j )
        coeff*= 2.0/n
        coeffs.append( coeff )
    print( "Chebychev coeffs", coeffs )

    print( "xnode, function at xnode, approximation" )
    for k in range( n ):
        xnode= xnodes[k]
        tnode= tnodes[k]
        p= 0.0
        for j in range( n ):
            p+= coeffs[j]*chbchvpoly( tnode, j )
        p-= 0.5*coeffs[0]
        print( xnode, fun( xnode ), p )
     
    return

# Interpolate Gauss with fitted Chebychev from RooFit
def gauss_test( n=5, nevent=10000 ):

    # Range and bins
    x= RooRealVar( "x", "x", -10, 10 )
    x.setBins( 25 )

    # Gauss properties
    mean= RooRealVar( "mean", "mean", 0, -10.0, 10.0 )
    sigma= RooRealVar( "sigma", "sigma", 5, 0.1, 10.0 )
    gauss= RooGaussian( "gauss", "signal p.d.f.", x, mean, sigma )

    # Generate a dataset from Gauss
    dataset= gauss.generate( RooArgSet( x ), nevent )
    datahist= dataset.binnedClone() 
    histpdf= RooHistPdf( "histpdf", "histpdf", RooArgSet( x ), datahist, 2 ) 

    # Return values from histogram
    def func( xx ):
        x.setVal( xx )
        return histpdf.getVal()

    # Setup Chebychev interpolation with python class above
    ca= ChebychevApproximation( func, -10.0, 10.0, n )
    #ca= ChebychevApproximation( func, -9.0, 9.0, n )
    coeffs= ca.getCoefficients()
    roocoeffs= [ coeff/coeffs[0] for coeff in coeffs ]
    print( roocoeffs )

    # Do the comparisons in RooFit, pure interpolation first, fitted second
    ral1= RooArgList( "ral1" )
    rrvs1= []
    for i in range( 1, len(coeffs) ):
        name= "1c"+str(i)
        title= "chbchv1 coeff "+str(i)
        rrv1= RooRealVar( name, title, roocoeffs[i], -10000.0, 10000.0 )
        rrv1.Print()
        ral1.add( rrv1 )
        rrvs1.append( rrv1 )        
    ral1.Print()
    chbchv1= RooChebychev( "chbchv1", "Chebychev p.d.f.", x, ral1 )
    chbchv1.Print()
    ral2= RooArgList( "ral2" )
    rrvs2= []
    for i in range( 1, len(coeffs) ):
        name= "2c"+str(i)
        title= "chbchv2 coeff "+str(i)
        rrv2= RooRealVar( name, title, roocoeffs[i], -10000.0, 10000.0 )
        rrv2.Print()
        ral2.add( rrv2 )
        rrvs2.append( rrv2 )        
    ral2.Print()
    chbchv2= RooChebychev( "chbchv2", "Chebychev p.d.f.", x, ral2 )
    chbchv2.Print()

    # Calculate chi^2 w/o or w/ fit
    cmdList= RooLinkedList( 1 )
    rca= RooCmdArg( RooFit.PrintLevel(0) )
    cmdList.Add( rca )
    chi2approx= chbchv1.createChi2( datahist ).getVal()        
    chbchv2.chi2FitTo( datahist, cmdList )
    chi2fit= chbchv2.createChi2( datahist ).getVal()
    print( "Chi^2 Chebychev approx vs data histogram", chi2approx )
    print( "Chi^2 Chebychev fit vs data histogram", chi2fit )

    # RooFit plots
    xframe= x.frame( RooFit.Title( "Gauss model Chebychev interpolation") )
    dataset.plotOn( xframe )
    chbchv1.plotOn( xframe, LineColor( kBlue ), Name( "chbchv1" ) )
    chbchv2.plotOn( xframe, LineColor( kRed ), Name( "chbchv2" ) )
    gauss.plotOn( xframe, LineColor( kGreen ), Name( "gauss" ) )
    xframe.Draw()
    global legend    
    legend= TLegend( 0.3, 0.2, 0.6, 0.4 )
    legend.AddEntry( xframe.findObject( "gauss" ), "Gauss(0,5)", "l" )
    legend.AddEntry( xframe.findObject( "chbchv1" ), "Chebychev before fit", "l" )
    legend.AddEntry( xframe.findObject( "chbchv2" ), "Chebychev after fit", "l" )
    legend.Draw( "same" )
    
    return


# For a RooDataHist, do a Chebychev parametrisation and
# then fit to the data to optimise the coefficients
def parametriseDataHist( datahist, x, n=9, nip=3, opt="f" ):

    # RooHistPdf for interpolated pdf values
    histpdf= RooHistPdf( "histpdf", "histpdf", RooArgSet( x ), datahist, nip )

    # Calculate Chebychev coefficients and create RooChebychev pdf
    chbchv= parametrisePdf( histpdf, x, n )
    
    # Fit coefficients if option set
    if "f" in opt:
        cmdList= RooLinkedList( 1 )
        rca= RooCmdArg( RooFit.PrintLevel(1) )
        cmdList.Add( rca )
        chbchv.chi2FitTo( datahist, cmdList )
    chi2= chbchv.createChi2( datahist ).getVal()
    print( "parametriseDataHist: Chi^2 Chebychev vs data histogram", chi2 )

    # The End
    return chbchv


# Calculate the Chebychev parametrisation of a RooPdf
# (or anything else with a working getVal()) with a
# RooChebychev (these work with normalised Chebychev coeficients)
def parametrisePdf( pdf, var, n=9 ):

    # Function for ChebychevApproximation
    def func( xx ):
        var.setVal( xx )
        return pdf.getVal()

    # Initialise approximation and get coefficients, convert to normalised
    # RooChebychev coefficients and create RooChebychev instance
    varrange= var.getRange()
    a= varrange.first
    b= varrange.second
    ca= ChebychevApproximation( func, a, b, n )
    coeffs= ca.getCoefficients()
    roocoeffs= [ coeff/coeffs[0] for coeff in coeffs ]
    rooCoeffsRal= RooArgList( "rcral" )
    for icoeff in range( 1, len(coeffs) ):
        name= "c"+str(icoeff)
        title= "chbchv coeff "+str(icoeff)
        roocoeff= roocoeffs[icoeff]
        rooCoeffRrv= RooRealVar( name, title, roocoeff )
        rooCoeffRrv.setError( abs(roocoeff)/3.0 )
        rooCoeffRrv.setConstant( False )
        SetOwnership( rooCoeffRrv, False )
        rooCoeffsRal.add( rooCoeffRrv )
    chbchv= RooChebychev( "chbchv", "Chebychev p.d.f.", var, rooCoeffsRal )
    chbchv.Print()

    # The End
    return chbchv

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
    return hists

# Parametrise with fit, plots
def parametrisationPlots():    
    mlbhistKeys= [ "h_mlb_171.0", "h_mlb_172.0", "h_mlb_172.5", "h_mlb_173.0", "h_mlb_174.0" ]
    mlbhists= getHistosFromFile( mlbhistKeys, "output_signal.root" )
    mtopDependence( mlbhists, n=8, a=40.0, b=148.0, txt="ParametrisationRooFit_ttbar"  )
    mtopDependence( mlbhists, n=8, a=40.0, b=160.0, txt="ParametrisationRooFit_ttbar"  )
    mtopDependence( mlbhists, n=9, a=40.0, b=148.0, txt="ParametrisationRooFit_ttbar" )
    mtopDependence( mlbhists, n=9, a=40.0, b=160.0, txt="ParametrisationRoofit_ttbar"  )
    mtopDependence( mlbhists, n=10, a=40.0, b=148.0, txt="ParametrisationRooFit_ttbar" )
    mtopDependence( mlbhists, n=10, a=40.0, b=160.0, txt="ParametrisationRoofit_ttbar"  )
    # mlbhistKeys= [ "tReco_mlb_minavg_171.0", "tReco_mlb_minavg_172.0",
    #                "tReco_mlb_minavg_172.5", "tReco_mlb_minavg_173.0", 
    #                "tReco_mlb_minavg_174.0" ]
    # mlbhists= getHistosFromFile( mlbhistKeys, "output_signal.root" )
    # mtopDependence( mlbhists, n=8, a=40.0, b=148.0, txt="Parametrisation_WWbb"  )
    # mtopDependence( mlbhists, n=8, a=40.0, b=160.0, txt="Parametrisation_WWbb"  )
    # mtopDependence( mlbhists, n=9, a=40.0, b=148.0, txt="Parametrisation_WWbb" )
    # mtopDependence( mlbhists, n=9, a=40.0, b=160.0, txt="Parametrisation_WWbb"  )
    return
def mtopDependence( mlbhists, n=9, a=40.0, b=148.0, txt="Parametrisation", opt="f" ):

    # Variable for RooFit
    x= RooRealVar( "mlb", "mlb", a, b )
    mlbhistKeys= mlbhists.keys()
    key= list(mlbhistKeys)[0]
    hist= mlbhists[key]
    binw= hist.GetBinWidth( 1 )
    nbins= int((b-a)/binw)
    print( "mtopDependence: bins in range", a, "to", b, "from", key, ":", nbins )
    x.setBins( nbins )
    
    # Do the parametrisation of the hists
    rooCoeffRrvsDict= dict()
    chbchvPdfs= dict()
    datahists= dict()
    for key in mlbhistKeys:
        mlbhist= mlbhists[key]
        datahist= RooDataHist( key, key, RooArgList( x ), mlbhist )
        chbchvPdf= parametriseDataHist( datahist, x, n=n, opt=opt )
        chbchvPdf.SetName( "chbchv_"+key )
        rooCoeffRrvs= chbchvPdf.getParameters( RooArgSet( x ) )
        rooCoeffRrvs.sort()
        rooCoeffRrvsDict[key]= rooCoeffRrvs
        chbchvPdfs[key]= chbchvPdf
        datahists[key]= datahist
    print( "mtopDependence: Results for RooChebychev coefficients" )
    for key in sorted( mlbhistKeys ):
        rooCoeffRrvs= rooCoeffRrvsDict[key]
        for rooCoeffRrv in rooCoeffRrvs:            
            value= rooCoeffRrv.getVal()
            error= rooCoeffRrv.getError()
            print( "{0:7.4f} {1:7.4f}".format( value, error ), end=" " )
        print()

    # Some globals so plots stay alive
    global canv1, canv2, canv3, tges, TH1Dchbchvhists, TH1Ddatahists, ratioPlots

    # Plots of parametrisations
    canv1= TCanvas( "canv1", "Fit plots", 600, 800 )
    canv1.Divide( 2, 3 )
    icanv= 0
    for key in mlbhistKeys:
        icanv+= 1
        canv1.cd( icanv )
        xframe= x.frame( RooFit.Title( "Mlb Chebychev "+key ) )
        datahist= datahists[key]
        datahist.plotOn( xframe )
        chbchv= chbchvPdfs[key]
        chbchv.plotOn( xframe )
        chi2= chbchv.createChi2( datahist ).getVal()
        ndof= nbins-n+1
        chi2ndof= chi2/float(ndof)
        chi2str= "{0:5.2f}".format( chi2ndof )
        chbchv.paramOn( xframe, Label( "Chi^2/ndof = "+chi2str ),
                        Format( "NE", AutoPrecision(1) ) )
        xframe.Draw()
    pdffilename= txt+"_"+str(n)+"_"+str(a)+"_"+str(b)+".pdf"
    canv1.Print( pdffilename+"(", "pdf" )
    
    # Fit info on TGraphError plots
    gStyle.SetOptFit( 1111 )

    # Ratio plots
    canv2= TCanvas( "canv2", "Significance plots", 600, 800 )
    canv2.Divide( 2, 3 )
    canv2.SetGrid()
    icanv= 0
    TH1Dchbchvhists= list()
    TH1Ddatahists= list()
    ratioPlots= list()
    for key in mlbhistKeys:
        icanv+= 1
        pad= canv2.cd( icanv )
        pad.SetGrid()
        chbchvPdf= chbchvPdfs[key]
        datahist= datahists[key]
        TH1Ddatahist= datahist.createHistogram( "TH1Ddatahist"+key, x )
        TH1Ddatahists.append( TH1Ddatahist )
        events= TH1Ddatahist.GetSumOfWeights()
        chbchvBinned= chbchvPdf.generateBinned( RooArgSet( x ), events,
                                                RooFit.ExpectedData() )
        TH1Dchbchvhist= chbchvBinned.createHistogram( "TH1Dchbchvhist"+key, x )
        TH1Dchbchvhists.append( TH1Dchbchvhist )
        ratioPlot= TGraphErrors()
        ratioPlot.SetName( "ratio_"+key )
        ratioPlot.SetTitle( "significance "+key )
        ratioPlots.append( ratioPlot )
        for ibin in range( 1, TH1Ddatahist.GetNbinsX()+1 ):
            dataValue= TH1Ddatahist.GetBinContent( ibin )
            dataError= TH1Ddatahist.GetBinError( ibin )
            axisValue= TH1Ddatahist.GetBinCenter( ibin )
            fitValue= TH1Dchbchvhist.GetBinContent( ibin )
            ratio= (fitValue-dataValue)/dataError
            ratioPlot.SetPoint( ibin-1, axisValue, ratio )
            ratioPlot.SetPointError( ibin-1, 0.0, 1.0 )
        ratioPlot.SetMinimum( -3.0 )
        ratioPlot.SetMaximum( 4.0 )
        ratioPlot.GetXaxis().SetLimits( a, b )
        ratioPlot.SetMarkerStyle( 20 )
        ratioPlot.SetMarkerSize( 0.75 )
        ratioPlot.Fit( "pol0" )
        ratioPlot.Draw( "ap" )
    canv2.Print( pdffilename, "pdf" )

    # Plots of coefficients vs mtop
    tges= list()
    ncoeffs= len(rooCoeffRrvs)
    for icoeff in range( ncoeffs ):
        tge= TGraphErrors()
        tge.SetName( "tge"+str(icoeff+1) )
        tge.SetTitle( "Chebychev coeff c"+str(icoeff+1) )
        imtopPoint= 0
        for key in sorted( mlbhistKeys ):
            tokens= key.split("_")
            mtop= float( tokens[len(tokens)-1] )
            rooCoeffRrvs= rooCoeffRrvsDict[key]
            # No index into RooArgSet, but is sorted
            rooCoeffRRvList= [ coeff for coeff in rooCoeffRrvs ]
            rooCoeffRrv= rooCoeffRRvList[icoeff]
            value= rooCoeffRrv.getVal()
            error= rooCoeffRrv.getError()
            tge.SetPoint( imtopPoint, mtop, value )
            tge.SetPointError( imtopPoint, 0.0, error )
            imtopPoint+= 1
        tges.append( tge )
    canv3= TCanvas( "canv3", "Coefficients mtop dependence", 800, 800 )
    canv3.DivideSquare( n-1 )
    for icoeff in range( ncoeffs ):
        canv3.cd( icoeff+1 )
        tge= tges[icoeff]
        tge.SetMarkerColor( 1 )
        tge.SetMarkerStyle( 20 )
        tge.SetMarkerSize( 0.75 )
        tge.Fit( "pol1" )
        tge.Draw( "ap" )
    canv3.Print( pdffilename+")", "pdf" )
                
    return

# Parametrisation with plain PyROOT fits
def parametrisationPlots2():
    
    mlbhistKeys= [ "h_mlb_171.0", "h_mlb_172.0", "h_mlb_172.5", "h_mlb_173.0", "h_mlb_174.0" ]
    global mlbhists
    mlbhists= getHistosFromFile( mlbhistKeys, "output_signal.root" )
    mtopDependence2( mlbhists, n=8, a=40.0, b=148.0, txt="ParametrisationPy_ttbar"  )
    mtopDependence2( mlbhists, n=8, a=40.0, b=160.0, txt="ParametrisationPy_ttbar"  )
    mtopDependence2( mlbhists, n=9, a=40.0, b=148.0, txt="ParametrisationPy_ttbar" )
    mtopDependence2( mlbhists, n=9, a=40.0, b=160.0, txt="ParametrisationPy_ttbar"  )
    mtopDependence2( mlbhists, n=10, a=40.0, b=148.0, txt="ParametrisationPy_ttbar" )
    mtopDependence2( mlbhists, n=10, a=40.0, b=160.0, txt="ParametrisationPy_ttbar"  )
    
    # mlbhistKeys= [ "tReco_mlb_minavg_171.0", "tReco_mlb_minavg_172.0",
    #                "tReco_mlb_minavg_172.5", "tReco_mlb_minavg_173.0", 
    #                "tReco_mlb_minavg_174.0" ]
    # mtopDependence2( mlbhists, n=8, a=40.0, b=148.0, txt="Parametrisation2_WWbb"  )
    # mtopDependence2( mlbhists, n=8, a=40.0, b=160.0, txt="Parametrisation2_WWbb"  )
    # mtopDependence2( mlbhists, n=9, a=40.0, b=148.0, txt="Parametrisation2_WWbb" )
    # mtopDependence2( mlbhists, n=9, a=40.0, b=160.0, txt="Parametrisation2_WWbb"  )
    
    return

# Function object around histo for ChebychevApproximation
class HistFunc:
    def __init__( self, hist ):
        self.__hist= hist
    #def __call__( self, x ):
    #    return self.__hist.Interpolate( x )
    def __call__( self, x ):
        return self.__hist.GetBinContent( self.__hist.FindBin( x ) )
# Plain ROOT fit parametrisation
def parametriseHistogram( hist, a, b, n ):
    histFunc= HistFunc( hist )
    ca= ChebychevApproximation( histFunc, a, b, n )
    pars= ca.getCoefficients()
    histname= hist.GetName()
    fittf= TF1( histname+"_tf1", ca, a, b, n )
    chebtf= TF1( histname+"_cheb_tf1", ca, a, b, n )
    for i in range( n ):
        fittf.SetParameter( i, pars[i] )
        chebtf.SetParameter( i, pars[i] )
        parname= "c"+str(i)
        fittf.SetParName( i, parname )
        chebtf.SetParName( i, parname )
    fittf.Print()
    hist.Fit( fittf, "0", "", a, b )
    return fittf, chebtf

def mtopDependence2( mlbhists, n=9, a=40.0, b=160.0, txt="Parametrisation2", opt="f" ):

    # Prepare welcome printout
    mlbhistKeys= mlbhists.keys()
    key= list(mlbhistKeys)[0]
    hist= mlbhists[key]
    binw= hist.GetBinWidth( 1 )
    nbins= int((b-a)/binw)
    print( "mtopDependence: bins in range", a, "to", b, "from", key, ":", nbins )
    
    # Do the parametrisation of the hists
    global fittfs, chebtfs
    fittfs= dict()
    chebtfs= dict()
    for key in sorted( mlbhistKeys ):
        mlbhist= mlbhists[key]
        fittf, chebtf= parametriseHistogram( mlbhist, a, b, n )
        fittfs[key]= fittf
        chebtfs[key]= chebtf
    if "v" in opt:
        print( "mtopDependence2: Results for Chebychev coefficients" )
    hpars= dict()
    hparerrs= dict()
    for key in sorted( mlbhistKeys ):
        mlbhist= mlbhists[key]
        fittf= fittfs[key]
        pars= fittf.GetParameters()
        errs= fittf.GetParErrors()
        hpars[key]= pars
        hparerrs[key]= errs
        if "v" in opt:
            for i in range( n ):
                value= pars[i]
                error= errs[i]
                print( "{0:7.4f} {1:7.4f}".format( value, error ), end=" " )
            print()
        
    # Some globals so plots stay alive
    global canv, tges, ratioPlots

    # Fit info on plots
    gStyle.SetOptFit( 1111 )
    from ROOT import kBlue
    
    # Plots of parametrisations
    canv1= TCanvas( "canv1", "Fit plots", 600, 800 )
    canv1.Divide( 2, 3 )
    icanv= 0
    for key in sorted( mlbhistKeys ):
        icanv+= 1
        canv1.cd( icanv )
        mlbhist= mlbhists[key]
        mlbhist.Draw()
        chebtf= chebtfs[key]
        chebtf.SetLineColor( kBlue )
        chebtf.Draw( "same" )
        fittf= fittfs[key]
        fittf.Draw( "same" )
    pdffilename= txt+"_"+str(n)+"_"+str(a)+"_"+str(b)+".pdf"
    canv1.Print( pdffilename+"(", "pdf" )

    # Ratio plots
    canv2= TCanvas( "canv2", "Significance plots", 600, 800 )
    canv2.Divide( 2, 3 )
    icanv= 0
    ratioPlots= list()
    for key in sorted( mlbhistKeys ):
        icanv+= 1
        pad= canv2.cd( icanv )
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
    canv2.Print( pdffilename, "pdf" )
    
    # Plots of coefficients vs mtop
    tges= list()
    ncoeffs= n
    for icoeff in range( ncoeffs ):
        tge= TGraphErrors()
        tge.SetName( "tge"+str(icoeff+1) )
        tge.SetTitle( "Chebychev coeff c"+str(icoeff) )
        imtopPoint= 0
        for key in sorted( mlbhistKeys ):
            tokens= key.split("_")
            mtop= float( tokens[len(tokens)-1] )
            value= hpars[key][icoeff]
            error= hparerrs[key][icoeff]
            tge.SetPoint( imtopPoint, mtop, value )
            tge.SetPointError( imtopPoint, 0.0, error )
            imtopPoint+= 1
        tges.append( tge )
    canv3= TCanvas( "canv3", "Coefficients mtop dependence", 800, 800 )
    canv3.DivideSquare( n )
    for icoeff in range( ncoeffs ):
        canv3.cd( icoeff+1 )
        tge= tges[icoeff]
        tge.SetMarkerColor( 1 )
        tge.SetMarkerStyle( 20 )
        tge.SetMarkerSize( 0.75 )
        tge.Fit( "pol1" )
        tge.Draw( "ap" )
    canv3.Print( pdffilename+")", "pdf" )
                
    return



# Run as main
def main():
    parametrisationPlots2()
    return
if __name__ == '__main__':
   main()
