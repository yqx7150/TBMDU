%% The Code is created based on the method described in the following paper: 
%[1] Qiegen Liu, Kun Yang, Jianhua Luo, Yuemin Zhu, Dong Liang, Highly undersampled MRI reconstruction using two-level Bregmanized method 
%    with dictionary learning, IEEE Transactions on Medical Imaging, 2013, 32 (7): 1290-1301.
% Date : 06/7/2013 
% Version : 1.0 
% The code and the algorithm are for non-comercial use only. 
% Copyright 2013, Department of Electronic Information Engineering, Nanchang University. 
% The current version is not optimized. 


clear all; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
getd = @(p)path(path,p);% Add some directories to the path
getd('utils\');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                               %#######%%%%% step1 :generate data %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%#######%%%%% read sampling %%%%
load Q1.mat; mask = Q1;   
figure(356); imshow(fftshift(mask),[]);   %
n = size(mask,2);
fprintf(1, 'n=%d, k=%d, Unsamped=%f\n', n, sum(sum(mask)),1-sum(sum(mask))/n/n); %

%#######%%%%% read image %%%%
M0 = im2double(imread('t2axialbrain.jpg'));   
if (length(size(M0))>2);  M0 = rgb2gray(M0);   end      
if (max(M0(:))<2);   M0 = M0*255;    end
figure(456); imshow(M0,[]);        

%#######%%%%% generate K-data %%%%
sigma = 0;
y_noisy = mask.*(fft2(M0) + randn(n)*sigma*255 + (0+1i)*randn(n)*sigma*255); %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                               %#######%%%%% step2 :Reconstruction %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% parameter setting %%%%%%%%%%%%%%
lambda = 12; %
mu = 10;   %
dic_stepsize = 0.01; %step size of dictionary updating
iter = 25; %iter number
slidingDis = 2;   %overlapping patch sliding
bb=8; % block size     
RR=4; % redundancy factor

%%%%%%%%%%%%%%%%%%%% initialization %%%%%%%%%%%%%%
Im = real(ifft2(y_noisy));   % zero-filled result
figure(33);imagesc(Im);colormap(gray);axis off; axis equal;
Im_extensin = padarray(Im,[bb-1 bb-1],'symmetric');       %
[blkMatrix,idx] = my_im2col(Im_extensin,[bb,bb],slidingDis);
vecOfMeans = mean(blkMatrix);
blkMatrix  = blkMatrix-ones(size(blkMatrix,1),1)*vecOfMeans;
A = init_dictionary(bb,RR);   % init dictionary
Anorm      = get_operator_norm_Liu(A,0);   %
rho_without_alpha = 1./(Anorm*Anorm);
[atomsize, samplenum] = size(blkMatrix);
[atomsize, atomnum] = size(A);
X = zeros(atomnum,samplenum);
cm = zeros(atomsize, samplenum);
f = y_noisy;
temp_FFT = zeros(size(Im));
psnr_idex = zeros(1,iter);

%%%%%%%%%%%%%%%%%%%  main iteration£ºTwo-level iteration  %%%%%%%%%%%%%%%%%%%%%
% pls see page 6, Algorithm2 
for iteri = 1:iter  %outer iter
    disp(sprintf('iteri: %d ', iteri));  %
    beta = 1/180;  %
    for iterii = 1:5  %inner iter
        %%%%%%%%%%%%% step 01£º update image patch variables %%%%%%%%%%%%%%%
        for mm= 1:4 
            Y = (lambda*beta/(lambda+beta))*(blkMatrix-A*sparse(X)+(1/beta)*cm);  %%% Eq.(18)
            X = X + (rho_without_alpha/beta) * A'*Y;   %%% compute X
            X = soft(X,rho_without_alpha/beta);  %Eq. (20)
        end
        cm = Y ;        %%% compute Y
        A = A + dic_stepsize*Y*X';          %%% dictionary update, Eq. (16)
        %A = A ./ repmat( sqrt(sum(A.^2,1)), atomsize,1 );
        A = A ./ repmat( sqrt(sum(A.*conj(A),1)), atomsize,1 );  
        Anorm = get_operator_norm_Liu(A,0); 
        rho_without_alpha = 1./(Anorm*Anorm); %used in Eq.(20)  
        
        %%%%%%%%%%%%% step 02£º update image variables %%%%%%%%%%%%%%%
        blocks = (A*X +(1/lambda-1/beta)*Y)+ones(size(blkMatrix,1),1)*vecOfMeans;
        [IMout,Weight] = patch_RT(blocks,n+2*(bb-1),bb,idx);
        %%%% update in frequency domain, Eq.(26)
        temp_Im = IMout./Weight;  temp_Im = temp_Im(bb:n+(bb-1),bb:n+(bb-1));  %
        temp_Im(temp_Im>255)=255;  %if the intensity range is known
        temp_part = fft2(temp_Im);
        temp_FFT(mask==0) = temp_part(mask==0);  %
        temp_FFT(mask==1) = (beta*temp_part(mask==1) + mu*f(mask==1))/(beta+mu); %sampled data
        Im = real(ifft2(temp_FFT));   %back to image domain
        Im(Im<0)=0;Im(Im>255)=255;  %if the intensity range is known
        
        %%%%%%%%%%%%% step 03£º update the patch %%%%%%%%%%%%%%%
        Im_extensin = padarray(Im,[bb-1 bb-1],'symmetric');       %extension
        [blkMatrix,idx] = my_im2col(Im_extensin,[bb,bb],slidingDis);
        vecOfMeans = mean(blkMatrix);
        blkMatrix = blkMatrix-ones(size(blkMatrix,1),1)*vecOfMeans;
        beta = beta*1.05;     
    end
    figure(197);II = displayDictionary_nonsquare2(A,0);  %display the dictionary 
    %%%%%%%%%%%%% step 11£º update outer constraint %%%%%%%%%%%%%%%
    f = f + (y_noisy - mask.*(fft2(Im)));    %%% Eq.(11)
    
    %%%%%%%%%%%%% step 12£º display outer result %%%%%%%%%%%%%%%
    figure(33);imagesc(Im);colormap(gray);axis off; axis equal;  %Rec image
    %figure(333);imagesc(abs(Im/255-M0/255));axis off; axis equal;colorbar;
    highfritererror(iteri)=norm(imfilter(abs(Im)/max(abs(Im(:))),fspecial('log',15,1.5)) - imfilter(M0/max(M0(:)),fspecial('log',15,1.5)),'fro');
    psnr_idex(iteri) = psnr(Im,M0,255);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                               %#######%%%%% step3 :display and save %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(173); plot(psnr_idex,'rv-');
param2.Im = Im;
param2.PSNR = psnr_idex;
param2.HFEN = highfritererror;
param2.Dictionary = A;

