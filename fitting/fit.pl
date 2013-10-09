#!/opt/star/bin/perl -w
#
use strict;
#naming the directories: 
my $workDir      = "/Users/nicolas/Project/ups2013/code/";
# after $workdir, the structure is
# pdfOutput                        |txtOutput
# date (I use /MMDD)
# vP_0p01       |    vP_0p05       |vP_0p01     |   vP_0p05
#pt_3 | pt_3p5 | pt_4 | pt_4p5 (in each vP subdirectory)

my $pdfPrefix = "pdfOutput/2009/";
my $txtPrefix = "txtOutput/2009/";
my $outFigDirStep = $workDir.$pdfPrefix;
my $outDatDirStep = $workDir.$txtPrefix;
my $outFigDir;
my $outDatDir;		

if(! -e $outFigDirStep){ system("mkdir $outFigDirStep");}
if(! -e $outDatDirStep){ system("mkdir $outDatDirStep");}
print "Welcome to SuperFitter! ";
###

#------------------------ OUTPUT DIRECTORIES per samples
my @sampleIndex  = ("","1","2","3","4","5","6","7"); #Input data sample. 1: pp@7TeV-d0 data ; 2: pp@7TeV-d3 data; 3: PbPb@276 regit 4: pp reco TT 5: pp reco GG 6: zhen tree, GG 7: pp@2.76TeV data
my $samstart = 3;
my $samend   = 3;
my @choseSamples = ("","pp7tev-d0", "pp7tev-d3","pbpbRegit_cm","pbpbPPrecoTT", "pbpbPPrecoGG","pbpbZhen2011","pp2p76tev");
# ------------------------------------------------------
my $doNarrowMass = 1; 
my $choseFitParams = 0; #what you want 0: (1s, 2s, 3s) 1: (1s, 2s/1s; 3s/1s); 2: (1S, (2s+3s)/1s); 3:(1S, 2S/1S, 3S/1S, 3S/2S)
my @choseWhat2Fit = ("yield_m7","ratio1","ratio23","ratio32");
my $prefix          = $choseWhat2Fit[$choseFitParams];

my $ptMuStart = 3;#single muon pt > 4GeV/c cut, that no-one wants to change... that makes me sad.
my $ptMuEnd = 3;
my @chooseptMu = ("0","3","3.5","4","4.5");
my @outFigPrefix_ptMu = ("pt_0","pt_3","pt_3p5","pt_4","pt_4p5");
my $ptMuPrefix;

my $vProbStart = 1;
my $vProbEnd = 1; #1
my @choosevProb = ("","0.01","0.05");
my @outFigPrefix_vProb = ("","vP_0p01","vP_0p05");
my $vProbPrefix;
# chose if fix sigma1:
my @sigma1  = ("0","1"); # fix or not sigma1
my $sigstart = 0;
my $sigend   = 1;

my @useRef   = ("0","1","2","3"); # 0 free; 1: fix to MC 2: fixed to MB  NOOOOO !!! look in the code for doc.
my $refstart  = 3;
my $refend    = 3;

# chose FSR param
#0: free; 1: both fixed 2: alpha fixed 3: npow fixed 
# keep in mind the shape parameters might not be the same for all bins... stick to the usual n=2.3 for convenience, for now.
my @fsr     = ("0","1","2","3");
my $fsrstart = 3;
my $fsrend   = 3;

my $centstart =0 ;
my $centend   =0 ;
#my @centrality_min = ("0","4", "8", "0", "50", "20","0");
#my @centrality_max = ("4","8", "24", "8","100", "40","0");
# ----------------------0----1---2---3----4----5----6----7---8---9---10---11--12---13--; 0->7 1S binning, 8->11 2S binning, n>11 for trials
 my @centrality_min = ("0", "0","2","4", "8","12","16","20","0","4", "8","0","16","28"); 
 my @centrality_max = ("40","2","4","8","12","16","20","28","4","8","40","8","28","40");
# 0-5 | 5-10 | 10-20 | 20-30 | 30-40 | 40-50(60) | 50(60)-100 | 0-20 | 20-90(alice comparisons)

my $modelstart =1;
my $modelend   =6;
#Background Model.  1: LS erf*exp + pol2; 2: LS RookeyPdf + pol2; 3: erf*exp; 4: pol2; 5: erf*exp+pol2 6: poly 3
my @bkgModel        = ("","1","2","3","4","5","6");

# fitting settings and plotting
my $paramsOnFigure = 1; #1: plot parameters;   0: plot CMS label
my $plotBkg        = 1; #0: hide LS or trkRot; 1: plot LS or trkRot data points and fit lines;
my $doTrkRot       = 0; #0: use LS;   1: use track rotation
my $doConstrainFit = 0; # 1: use constrain method

my $muonEtaMin     = -2.4;
my $muonEtaMax     = 2.4; 
my $dimuYMin       = -2.4; 
my $dimuYMax       = 2.4; #used when not binning in y_ups
my $upsPtCutMin    = 0;
my $upsPtCutMax    = 150; #used when not binning in pT
my $isam;

# my @rapBinMin = ("-2.4","-1.6","-0.8","0","0.8","1.6");
# my @rapBinMax = ("-1.6","-0.8","0","0.8","1.6","2.4");

# ----------------0----1-------2------3----4-----5-----6-----7--# rapidity binning suited for PbPb 1S (0->4), 2S (5->7);
my @rapBinMin = ("0.","0.4" ,"0.7","1.0","1.5","0." ,"1.2","0.");
my @rapBinMax = ("0.4","0.7","1.0","1.5","2.4","1.2","2.4","2.4");


# pT binning for 2S:
#my @upsPtBinMin = ("0","6.5","10","0");	 
#my @upsPtBinMax = ("6.5","10","20","20");
# ------------------0-----1-----2---3----4----5-----6-----7------8---9## pT binning for 1S (0->5), 2S (6->9);
my @upsPtBinMin = ("0"  ,"2.5","5","8" ,"12","20" ,"0"  ,"6.5","10","0");	 
my @upsPtBinMax = ("2.5","5" , "8","12","20","150","6.5","10","20","20");
my $dontDoRapNow =0;
my $dontDoPtNow =0;
#loops for mkdir purposes
for (my $ivProb=$vProbStart; $ivProb<=$vProbEnd; $ivProb++)
{
    $vProbPrefix = $outFigPrefix_vProb[$ivProb];
    my $outFigDirProb = $outFigDirStep.$vProbPrefix."/";
    my $outDatDirProb = $outDatDirStep.$vProbPrefix."/";
    
    if(! -e $outFigDirProb){ system("mkdir $outFigDirProb");}
    if(! -e $outDatDirProb){ system("mkdir $outDatDirProb");}
    
    for (my $iPtMu=$ptMuStart; $iPtMu<=$ptMuEnd; $iPtMu++)
    {
	$ptMuPrefix = $outFigPrefix_ptMu[$iPtMu];
	$outFigDir = $outFigDirProb.$ptMuPrefix."/";
	$outDatDir = $outDatDirProb.$ptMuPrefix."/";

	if(! -e $outFigDir){ system("mkdir $outFigDir");}
	if(! -e $outDatDir){ system("mkdir $outDatDir");}
    }
}

for( my $icent=$centstart; $icent <=$centend; $icent++)  
  {
    for ( $isam=$samstart; $isam<=$samend; $isam++)
      {
	$doNarrowMass = 0; 
	if ($isam==1) {$doNarrowMass=1;}
	#loop over vertex probabilities
	for (my $ivProb=$vProbStart; $ivProb<=$vProbEnd; $ivProb++)
	{
	    $vProbPrefix = $outFigPrefix_vProb[$ivProb];
	    my $vProb = $choosevProb[$ivProb];
	    # these two lines are repeated because $outFigDir and $outDatDir go in the function's arguments
	    my $outFigDirProb = $outFigDirStep.$vProbPrefix."/";
	    my $outDatDirProb = $outDatDirStep.$vProbPrefix."/";
    
	    for (my $iPtMu=$ptMuStart; $iPtMu<=$ptMuEnd; $iPtMu++)
	    {
		my $muonPtCut      = $chooseptMu[$iPtMu]; #single muon pT cut
     	        $ptMuPrefix = $outFigPrefix_ptMu[$iPtMu];	
		# these two lines are repeated because $outFigDir and $outDatDir go in the function's arguments
		$outFigDir = $outFigDirProb.$ptMuPrefix."/";
		$outDatDir = $outDatDirProb.$ptMuPrefix."/";
		#loop over fixed|unfixed parameters for systematics
		for (my $isig=$sigstart; $isig<=$sigend; $isig++)
		{
		    for (my $ifsr=$fsrstart; $ifsr<=$fsrend; $ifsr++)
		    {
			for (my $iref=$refstart; $iref<=$refend; $iref++)
			{
			    # for systm. studies
			    next if($iref==0 && ($isig!=0 && $ifsr!=0));
			    #  next if($iref==1);# do not fix to MC 
			    next if(($iref==1 || $iref==2) && ($isig==0 && $ifsr==0));
			    
			    for ( my $ibkg=$modelstart; $ibkg<=$modelend; $ibkg++)
			    {
				next if(($dontDoRapNow==1 && $dontDoPtNow==0)||($dontDoRapNow==0 && $dontDoPtNow==1) || (($dontDoRapNow==0 && $dontDoPtNow==0)));
				my $centMin = $centrality_min[$icent];
				my $centMax = $centrality_max[$icent];
				my $myfsr   = $fsr[$ifsr];
				my $myref   = $useRef[$iref];
				my $bkgType = $bkgModel[$ibkg];
				
				my $mysigma1= $sigma1[$isig];
				
				if ($ibkg == 1 || $ibkg==2) {$plotBkg=1;}
				else {$plotBkg=0;}
				
				my $choseSample = $sampleIndex[$isam];
				my $sample = $choseSamples[$choseSample];
				$upsPtCutMin=0;
				$upsPtCutMax=150;
				$dimuYMin = -2.4;
				$dimuYMax = 2.4;
				
				print"Cent: $centMin - $centMax \n BkgModel: $bkgType \n Sample: $sample \n Muon cut:	$muonEtaMin - $muonEtaMax\n";
				print "fsr:$myfsr \n sig1: $mysigma1 \n reference: $myref\n";
				print " Sample is: $sample\n";
				print  "vProb is: $vProb\n";
				print   "single Muon pT > $muonPtCut\n";
				print "$upsPtCutMin > pT_Ups > $upsPtCutMax";
				print "$dimuYMin > y_Ups > $dimuYMax";
				system("root -l -b -q	'fitUpsilonYields.C($choseSample,$choseFitParams,$bkgType,$myfsr,$mysigma1,$centMin,$centMax,$muonEtaMin,$muonEtaMax,$dimuYMin,$dimuYMax,$muonPtCut,$plotBkg,$doTrkRot,$doConstrainFit,$myref,$paramsOnFigure,\"$sample\",\"$outFigDir\",\"$outDatDir\",\"$prefix\",$doNarrowMass,$vProb)'
			"); 
			    }	#bkg 
			} # ref
		    }#fsr
		}# sigma
	    }
	}
      }#sample
  }#centrality loop
for ( $isam=$samstart; $isam<=$samend; $isam++)
{
    $doNarrowMass = 0; 
    if ($isam==1) {$doNarrowMass=1;}
    #loop over vertex probabilities
    for (my $ivProb=$vProbStart; $ivProb<=$vProbEnd; $ivProb++)
    {
	$vProbPrefix = $outFigPrefix_vProb[$ivProb];
	my $vProb = $choosevProb[$ivProb];
	# these two lines are repeated because $outFigDir and $outDatDir go in the function's arguments
	my $outFigDirProb = $outFigDirStep.$vProbPrefix."/";
	my $outDatDirProb = $outDatDirStep.$vProbPrefix."/";
	
	for (my $iPtMu=$ptMuStart; $iPtMu<=$ptMuEnd; $iPtMu++)
	{
	    my $muonPtCut      = $chooseptMu[$iPtMu]; #single muon pT cut
	    $ptMuPrefix = $outFigPrefix_ptMu[$iPtMu];	
	    # these two lines are repeated because $outFigDir and $outDatDir go in the function's arguments
	    $outFigDir = $outFigDirProb.$ptMuPrefix."/";
	    $outDatDir = $outDatDirProb.$ptMuPrefix."/";
	    #loop over fixed|unfixed parameters for systematics
	    for (my $isig=$sigstart; $isig<=$sigend; $isig++)
	    {
		for (my $ifsr=$fsrstart; $ifsr<=$fsrend; $ifsr++)
		{
		    for (my $iref=$refstart; $iref<=$refend; $iref++)
		    {
			# for systm. studies
			next if($iref==0 && ($isig!=0 && $ifsr!=0));
			#  next if($iref==1);# do not fix to MC 
			next if(($iref==1 || $iref==2) && ($isig==0 && $ifsr==0));
			
			for(my $iRap=0; $iRap<=7; $iRap++){ #hardcoded... so what?
			    $dimuYMin=$rapBinMin[$iRap];
			    $dimuYMax=$rapBinMax[$iRap];
			    $upsPtCutMin    = 0;
			    $upsPtCutMax    = 150; #used when not binning in pT
			    next if($dontDoRapNow==1);
			    for ( my $ibkg=$modelstart; $ibkg<=$modelend; $ibkg++)
			    {
				
				my $centMin = $centrality_min[0];
				my $centMax = $centrality_max[0];
				my $myfsr   = $fsr[$ifsr];
				my $myref   = $useRef[$iref];
				my $bkgType = $bkgModel[$ibkg];
				
				my $mysigma1= $sigma1[$isig];
				
				if ($ibkg == 1 || $ibkg==2) {$plotBkg=1;}
				else {$plotBkg=0;}
				
				my $choseSample = $sampleIndex[$isam];
				my $sample = $choseSamples[$choseSample];
				
				
				print"Cent: $centMin - $centMax \n BkgModel: $bkgType \n Sample: $sample \n Muon cut:	$muonEtaMin - $muonEtaMax\n";
				print "fsr:$myfsr \n sig1: $mysigma1 \n reference: $myref\n";
				print " Sample is: $sample\n";
				print  "vProb is: $vProb\n";
				print   "single Muon pT > $muonPtCut\n";
				print    " $upsPtCutMin < pT_Ups < $upsPtCutMax \n";
				print     "  $dimuYMin < y_Ups < $dimuYMax \n";
				system("root -l -b -q	'fitUpsilonYields_Yvariant.C($choseSample,$choseFitParams,$bkgType,$myfsr,$mysigma1,$centMin,$centMax,$muonEtaMin,$muonEtaMax,$dimuYMin,$dimuYMax,$muonPtCut,$plotBkg,$doTrkRot,$doConstrainFit,$myref,$paramsOnFigure,\"$sample\",\"$outFigDir\",\"$outDatDir\",\"$prefix\",$doNarrowMass,$vProb)'
			"); 
			    } #bkg 
			} # rapidity bins
			for(my $iUpsPt=0; $iUpsPt<=9; $iUpsPt++){ #hardcoded... so what?
			    $upsPtCutMin=$upsPtBinMin[$iUpsPt];
			    $upsPtCutMax=$upsPtBinMax[$iUpsPt];
			    $dimuYMin       = -2.4; 
			    $dimuYMax       = 2.4; #used when not binning in y_ups
			    next if($dontDoPtNow==1);
			    for ( my $ibkg=$modelstart; $ibkg<=$modelend; $ibkg++)
			    {
				
				my $centMin = $centrality_min[0];
				my $centMax = $centrality_max[0];
				my $myfsr   = $fsr[$ifsr];
				my $myref   = $useRef[$iref];
				my $bkgType = $bkgModel[$ibkg];
				
				my $mysigma1= $sigma1[$isig];
				
				if ($ibkg == 1 || $ibkg==2) {$plotBkg=1;}
				else {$plotBkg=0;}
				
				my $choseSample = $sampleIndex[$isam];
				my $sample = $choseSamples[$choseSample];
				
				
				print"Cent: $centMin - $centMax \n BkgModel: $bkgType \n Sample: $sample \n Muon cut:	$muonEtaMin - $muonEtaMax\n";
				print "fsr:$myfsr \n sig1: $mysigma1 \n reference: $myref\n";
				print " Sample is: $sample\n";
				print  "vProb is: $vProb\n";
				print   "single Muon pT > $muonPtCut\n";
				print "$upsPtCutMin > pT_Ups > $upsPtCutMax";
				print "$dimuYMin > y_Ups > $dimuYMax";
				system("root -l -b -q	'fitUpsilonYields_variant.C($choseSample,$choseFitParams,$bkgType,$myfsr,$mysigma1,$centMin,$centMax,$muonEtaMin,$muonEtaMax,$upsPtCutMin,$upsPtCutMax,$muonPtCut,$plotBkg,$doTrkRot,$doConstrainFit,$myref,$paramsOnFigure,\"$sample\",\"$outFigDir\",\"$outDatDir\",\"$prefix\",$doNarrowMass,$vProb)'
			"); 
			    } #bkg 
			} # pT bins
		    } # ref
		}#fsr
	    }# sigma
	}
    }
}#sample
