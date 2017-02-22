clear all
close all
clc



reorderBlocks = [4 1 3 6 5 2];
filename_Gal = 'CUST-6x7K-160602-addRiboswitch-5rep-JeffreyMunzar.csv';
filename_Hyb = 'AfterHyb_HighMg_g20xdr20_G_SLOT02_S01_L transposed rotated 180_30pixel.csv';
filename_Assay = 'AfterAssay_HighMg_g20xdr20_G_SLOT02_S01_L transposed rotated 180_30pixel.csv';
filename_ACEsequences = 'CUST-6x7K-160602-addRiboswitch-5rep-JeffreyMunzar-SequencesDeltaGs.csv';
Gain_Hyb = 4;
Gain_Assay = 4;
minACEsAfterOutlierRemoval = 3;
aptamerSequenceLength = 87;
misMatch = [10 1];
misMatchIndex = [554 1333 1334 1334]; 
startingTilingInAptamer = 1;
workingfolder = 'addRiboswitch_4_10C_12mM';
Flag_SkipKm = 0;
LigandConc = [.004 .0004 .00004 .000004]; % Adenine concentrations, Molar
WashingBlank1 = 0; 
WashingBlank2 = 100;
HybLim1 = 0; 
HybLim2 = 20000;
Baseline1 = 0; 
Baseline2 = 75;
Ligand1 = 0; 
Ligand2 = 25;
Score1_1 = 0.2;
Score1_2 = 5;
Score2_1 = 1;
Score2_2 = 8;
QuantitativeFlag = 0;
SurfaceDensityFlag = 0;
CompareRedGreenChannelsFlag = 0;
CompareSodiumPotassiumFlag = 0;

RiboswitchFlag = 1;
RiboswitchFlagCovariation = 0

filename_CovationData = 'covariation_ImportToMatlab.csv'



