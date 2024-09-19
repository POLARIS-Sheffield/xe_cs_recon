function [OutTable] = xe_cs_recon(kdata, ktraj, mtx_reco, flipDP, flipGP, TE, T2TP, T2GAS, AF, saveStatus)

    %% Dissolved 129Xe Compressed Sensing Reconstruction  
    
    % ------------------------------------------------------------------------%
    
    % Inputs:
    
    % kdata --- [3, nsamples, nspokes], k-space data (Gas, M, RBC)
    % ktraj --- [3, nsamples, nspokes], k-space trajectory
    % mtx_reco --- reconstructed matrix size
    % flipDP --- flip angle on the dissolved-phase (22deg)
    % flipGP --- flip angle on the gas-phase (0.22deg)
    % TE --- echo time
    % T2TP --- T2* of the membrane signal, from calibration spectra
    % T2GAS -- T2* of the gas signal, from calibration spectra
    % AF -- acceleration factor
    % saveStatus --- 0 = do not save, 1 = save
    
    % ------------------------------------------------------------------------%
    
    % Outputs:
    
    % OutTable --- table containing SNR, mean, median, std, min, max and CV for RBC:M, RBC:Gas and M:Gas   
    
    % ------------------------------------------------------------------------%
    
    % Contributors:
    
    % Jemima H. Pilgrim-Morris
    % Guilhem J. Collier
    
    % ------------------------------------------------------------------------%

    % Literature:

    % https://github.com/mrirecon/bart-webinars/blob/master/webinar1/day1_advanced_recon.ipynb
    % Sparse MRI: The application of compressed sensing for rapid MR imaging, Lustig et al, MRM, 2007, https://doi.org/10.1002/mrm.21391

    %%
    
    if nargin < 10
        saveStatus = 0;
    end

    if saveStatus == 1
        mkdir('CS_recon')
        cd CS_recon
    end


    %% Undersample

    nspokes = size(ktraj, 3);

    X = round(nspokes/AF);
    
    ktraj = ktraj(:,:,1:X);
    kdata_GAS = kdata(1,:,1:X);
    kdata_M = kdata(2,:,1:X);
    kdata_RBC = kdata(3,:,1:X);


    %% M Recon
    
    % inverse NUFFT
    nufft_m1 = bart(['nufft -i -x', num2str(mtx_reco),':',num2str(mtx_reco),':',num2str(mtx_reco)], ktraj, kdata_M);
    
    % compute Cartesian k-space
    kspM = bart('fft -u 3', nufft_m1);
    
    % compute sensitivities using ESPIRiT
    sensM = bart('ecalib -m1', kspM);
    
    % CS recon
    CS_imgM1 = bart('pics -m -R I:0:0.003 -R T:7:0:0.0003 -S -t', ktraj, kdata_M, sensM);
    CS_imgM = bart('flip 3', CS_imgM1);
    M_CS_abs = double(abs(CS_imgM));
    
    % create mask
    N_M_CS =real(CS_imgM(:,:,[3:6 (size(CS_imgM,3)-5):(size(CS_imgM,3)-2)]));
    noise_M_CS =std(N_M_CS(:));
    Mask_CS = logical(zeros(mtx_reco, mtx_reco, mtx_reco));
    Mask_CS(M_CS_abs>(150*noise_M_CS))=1;
    Mask_CS=bwareaopen(Mask_CS,120,18);
    signal_M_CS = mean(M_CS_abs(Mask_CS));
    SNR_M_CS = signal_M_CS/noise_M_CS;
    
    figure
    montage(Mask_CS, Indices=10:22, ThumbnailSize=[]);
    title('CS Mask')
    set(gca,'FontName','Arial','FontSize',15);
    if saveStatus == 1
        print('CS_Mask','-r300','-dpng');
    end

    % Plot M image
    figure
    montage(M_CS_abs, Indices=10:22, ThumbnailSize=[]);
    clim([0, max(M_CS_abs, [], 'all')])
    title('M')
    set(gca,'FontName','Arial','FontSize',15);
    if saveStatus == 1
        print('CS_M_img','-r300','-dpng');
    end

%% RBC Recon
    
    % inverse NUFFT
    nufft_rbc1 = bart(['nufft -i -x', num2str(mtx_reco),':',num2str(mtx_reco),':',num2str(mtx_reco)], ktraj, kdata_RBC);
    
    % compute Cartesian k-space
    kspR = bart('fft -u 3', nufft_rbc1);
    
    % compute sensitivities using ESPIRiT
    sensR = bart('ecalib -m1', kspR);
    
    % CS recon
    CS_imgR1 = bart('pics -m -R I:0:0.003 -R T:7:0:0.0003 -S -t', ktraj, kdata_RBC, sensR);
    CS_imgR = bart('flip 3', CS_imgR1);
    RBC_CS_abs = double(abs(CS_imgR));
    
    % SNR
    N_RBC_CS =real(CS_imgR(:,:,[3:6 (size(CS_imgR,3)-5):(size(CS_imgR,3)-2)]));
    noise_RBC_CS =std(N_RBC_CS(:));
    signal_RBC_CS = mean(RBC_CS_abs(Mask_CS));
    SNR_RBC_CS = signal_RBC_CS/noise_RBC_CS;
    
    % Plot RBC image
    figure
    montage(RBC_CS_abs, Indices=10:22, ThumbnailSize=[]);
    clim([0, max(RBC_CS_abs, [], 'all')])
    title('RBC')
    set(gca,'FontName','Arial','FontSize',15);
    if saveStatus == 1
        print('CS_RBC_img','-r300','-dpng');
    end
    
    %% Gas Recon
    
    % inverse NUFFT
    nufft_g1 = bart(['nufft -i -x', num2str(mtx_reco),':',num2str(mtx_reco),':',num2str(mtx_reco)], ktraj, kdata_GAS);
    
    % compute Cartesian k-space
    kspG = bart('fft -u 3', nufft_g1);
    
    % compute sensitivities using ESPIRiT
    sensG = bart('ecalib -m1', kspG);
    
    % CS recon
    CS_imgG1 = bart('pics -m -R I:0:0.003 -R T:7:0:0.0003 -S -t', ktraj, kdata_GAS, sensG);
    CS_imgG = bart('flip 3', CS_imgG1);
    G_CS_abs = double(abs(CS_imgG))*sind(flipDP)/sind(flipGP)*exp(-TE/T2TP)/exp(-TE/T2GAS);
    
    % SNR
    N_G_CS =real(CS_imgG(:,:,[3:6 (size(CS_imgG,3)-5):(size(CS_imgG,3)-2)]));
    noise_G_CS =std(N_G_CS(:));
    signal_G_CS = mean(G_CS_abs(Mask_CS));
    SNR_G_CS = signal_G_CS/noise_G_CS;
    
    % Plot Gas image
    figure
    montage(G_CS_abs, Indices=10:22, ThumbnailSize=[]);
    clim([0, max(G_CS_abs, [], 'all')])
    title('Gas')
    set(gca,'FontName','Arial','FontSize',15);
    if saveStatus == 1
        print('CS_Gas_img','-r300','-dpng');
    end
    
    %% Ratio Maps

    % RBC:M
    RBC_M_CS = (RBC_CS_abs./M_CS_abs).*Mask_CS;
    
    figure
    montage(RBC_M_CS, Indices=10:22, ThumbnailSize=[]);
    clim([0, max(RBC_M_CS, [], 'all')])
    colorbar
    colormap("hot")
    title('RBC:M')
    set(gca,'FontName','Arial','FontSize',15);
    set(gcf, 'Position', get(0, 'Screensize'));
    if saveStatus == 1
        print('CS_RBCM_map','-r300','-dpng');
    end

    
    % RBC:GAS
    RBC_GAS_CS = (RBC_CS_abs./G_CS_abs).*Mask_CS;
    
    figure
    montage(RBC_GAS_CS, Indices=10:22, ThumbnailSize=[]);
    clim([0, max(RBC_GAS_CS, [], 'all')])
    colorbar
    colormap("hot")
    title('RBC:Gas')
    set(gca,'FontName','Arial','FontSize',15);
    set(gcf, 'Position', get(0, 'Screensize'));
    if saveStatus == 1
        print('CS_RBCGas_map','-r300','-dpng');
    end
    
    % M:GAS
    M_GAS_CS = (M_CS_abs./G_CS_abs).*Mask_CS;
    
    figure
    montage(M_GAS_CS, Indices=10:22, ThumbnailSize=[]);
    clim([0, max(M_GAS_CS, [], 'all')])
    colorbar
    colormap("hot")
    title('M:Gas')
    set(gca,'FontName','Arial','FontSize',15);
    set(gcf, 'Position', get(0, 'Screensize'));
    if saveStatus == 1
        print('CS_MGas_map','-r300','-dpng');
    end

    %% Distribution metrics
    
    % RBC:M
    meanRBCMCS = mean(RBC_M_CS(Mask_CS));
    medRBCMCS = median(RBC_M_CS(Mask_CS));
    stdRBCMCS = std(RBC_M_CS(Mask_CS));
    minRBCMCS = min(RBC_M_CS(Mask_CS));
    maxRBCMCS = max(RBC_M_CS(Mask_CS));
    cvRBCMCS = stdRBCMCS/meanRBCMCS;
    
    % RBC:GAS
    meanRBCGASCS = mean(RBC_GAS_CS(Mask_CS));
    medRBCGASCS = median(RBC_GAS_CS(Mask_CS));
    stdRBCGASCS = std(RBC_GAS_CS(Mask_CS));
    minRBCGASCS = min(RBC_GAS_CS(Mask_CS));
    maxRBCGASCS = max(RBC_GAS_CS(Mask_CS));
    cvRBCGASCS = stdRBCGASCS/meanRBCGASCS;
    
    % M:GAS
    meanMGASCS = mean(M_GAS_CS(Mask_CS));
    medMGASCS = median(M_GAS_CS(Mask_CS));
    stdMGASCS = std(M_GAS_CS(Mask_CS));
    minMGASCS = min(M_GAS_CS(Mask_CS));
    maxMGASCS = max(M_GAS_CS(Mask_CS));
    cvMGASCS = stdMGASCS/meanMGASCS;
    
    %% Out Table
    
    OutTable = table(AF, SNR_RBC_CS, SNR_M_CS, SNR_G_CS, ...
        meanRBCMCS, medRBCMCS, stdRBCMCS, minRBCMCS, maxRBCMCS, cvRBCMCS,...
        meanRBCGASCS, medRBCGASCS, stdRBCGASCS, minRBCGASCS, maxRBCGASCS, cvRBCGASCS, ...
        meanMGASCS, medMGASCS, stdMGASCS, minMGASCS, maxMGASCS, cvMGASCS);
    if saveStatus == 1
        save('CS_recon')
        writetable(OutTable, 'CS_recon_results.xlsx');
    end


end