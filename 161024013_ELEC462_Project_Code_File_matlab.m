%%  ELEC462 FUNDAMENTALS OF STATISTICAL SIGNAL PROCESSING
%     TERM PROJECT
%     ESTIMATION SIGNAL WITH THREE DIFFERENT ESTIMATORS
%     Maximum Likelihood Estimation 
%     Maximum A-Posteriori Estimation 
%     Linear MMSE estimation

%     STUDENT : %% 
%     ID : %%
%     Github: https://github.com/vatozZ

close all; clear all; clc;
load('SigRes.mat');


% b) Maximum Likelihood Estimation 
% x_ML = argmax(log(f(y|x)))


figure(1);
x_ML = inv(transpose(H)*H) * transpose(H) * y; % Maximum Likelihood function
plot(x_ML,'LineWidth',3,...
    'MarkerSize',10,...
    'Color','r');
title('ML Estimated and Ground Truth Signal');
xlabel('samples');
ylabel('magnitude');
hold on;
plot(x,'LineWidth',3,...
    'MarkerSize',10,...
    'Color','b');
legend('x: maximum Likelihood estimation', 'x:  ground Truth signal');

% d) Maximum A-Posteriori Estimation 
% x_ML = argmax(log(f(y|x) . f(x)))

var_W = 25*(10^-4); % variance of Gaussian noise (W)
mean_W = 0; % noises mean
mean_X = 0; % X's mean 
cov_X = zeros(1,4); % covariance matrix of X will be fill in loop with changing 
                                % Lambda smoothing parameters
x_MAP = zeros(1,4); % X map results will be calculated
lambda_ = [1 10 100 1000] ; %lambda smoothing parameters
color_ =['r' 'y' 'b' 'g']; 
figure();
PSNR_dB_MAP_ = zeros(1, 4);
for i = 1:4
    if i == 1
        cov_X = 1 * (L.') * L ;    
        x_MAP = inv((H.') * H + var_W * cov_X) * (H.') *  y  +zeros(256);
        plot(x, 'Color', 'k','LineWidth',3, 'DisplayName', 'x ground truth');
        hold on ; 
        text_ = ['MAP for \lambda = ', num2str(lambda_(i))];
        title('MAP Estimated Signals via changing \lambda and Ground Truth Signal');
        plot(x_MAP, 'Color',color_(i),'LineWidth',2,'DisplayName',text_);
        xlabel('samples');
        ylabel('magnitude');
        hold on ;
        legend();
        PSNR_dB_MAP_1 = 20*log10(max(x)*sqrt(256)/(norm(x - x_MAP)));
    end
    if i == 2
        cov_X = 10 * (L.') * L ;    
        x_MAP = inv((H.') * H + var_W * cov_X) * (H.') *  y  +zeros(256);
        plot(x, 'Color', 'k','LineWidth',3, 'DisplayName', 'x ground truth');
        hold on ; 
        text_ = ['MAP for \lambda = ', num2str(lambda_(i))];
        title('MAP Estimated Signals via changing \lambda and Ground Truth Signal');
        plot(x_MAP, 'Color',color_(i),'LineWidth',2,'DisplayName',text_);
        xlabel('samples');
        ylabel('magnitude');
        hold on ;
        legend();
        PSNR_dB_MAP_10 = 20*log10(max(x)*sqrt(256)/(norm(x - x_MAP)));
    end
    if i == 3
        cov_X = 100 * (L.') * L ;    
        x_MAP = inv((H.') * H + var_W * cov_X) * (H.') *  y  +zeros(256);
        plot(x, 'Color', 'k','LineWidth',3, 'DisplayName', 'x ground truth');
        hold on ; 
        text_ = ['MAP for \lambda = ', num2str(lambda_(i))];
        title('MAP Estimated Signals via changing \lambda and Ground Truth Signal');
        plot(x_MAP, 'Color',color_(i),'LineWidth',2,'DisplayName',text_);
        xlabel('samples');
        ylabel('magnitude');
        hold on ;
        legend();
        PSNR_dB_MAP_100 = 20*log10(max(x)*sqrt(256)/(norm(x - x_MAP)));
    end
    if i == 4
        cov_X = 1000 * (L.') * L ;    
        x_MAP = inv((H.') * H + var_W * cov_X) * (H.') *  y  +zeros(256);
        plot(x, 'Color', 'k','LineWidth',3, 'DisplayName', 'x ground truth');
        hold on ; 
        text_ = ['MAP for \lambda = ', num2str(lambda_(i))];
        title('MAP Estimated Signals via changing \lambda and Ground Truth Signal');
        plot(x_MAP, 'Color',color_(i),'LineWidth',2,'DisplayName',text_);
        xlabel('samples');
        ylabel('magnitude');
        hold on ;
        legend();
        PSNR_dB_MAP_1000 = 20*log10(max(x)*sqrt(256)/(norm(x - x_MAP)));
    end

end
%x_LMMSE = cov_X*(H.')* (inv(H*cov_X*(H.') + var_W*var_W*identity_)) * ( y-H*mean_X)+ mean_X;

% e) For Linear MMSE estimation
% f(x|y)  = [f(y|x)f(x)]/[integral(f(y|x')f(x')dx'), (-inf,+inf)]

identity_ = eye(256); % I - identity ( 256x256)
 %already declared above, constant parameters
var_W = 25*(10^-4); 
mean_W = 0; 
mean_X = 0;
figure();

% same loop with above but calculated Linear Minimum Mean Square Error Estimator in this case.
for i = 1:4
    if i == 1
        cov_X = 1 * (L.') * L ;    
        x_LMMSE = cov_X * (H.')* (inv(H*cov_X*(H.') + var_W*var_W*identity_)) *y + zeros(256,1);
        plot(x, 'Color', 'k','LineWidth',4, 'DisplayName', 'x ground truth');
        hold on;
        plot(x_LMMSE, 'Color','r','LineWidth',1,'DisplayName','\lambda = 1');
        hold on;
        result_MMSE_1 = 20*log10(max(x)*sqrt(256)/(norm(x - x_LMMSE)));
    end
    
    if i== 2
        cov_X = 10 * (L.') * L ;    
        x_LMMSE = cov_X*(H.')* (inv(H*cov_X*(H.') + var_W*var_W*identity_)) *y + zeros(256,1);
        plot(x_LMMSE, 'Color','g','LineWidth',2,'DisplayName','\lambda = 10');
        hold on;
        result_MMSE_10 =20*log10(max(x)*sqrt(256)/(norm(x - x_LMMSE)));
    end
    
    if i== 3
        cov_X = 100 * (L.') * L ;    
        x_LMMSE = cov_X*(H.')* (inv(H*cov_X*(H.') + var_W*var_W*identity_)) * y +  zeros(256,1);
        plot(x_LMMSE, 'Color','b','LineWidth',2,'DisplayName','\lambda = 100');
        hold on;
        result_MMSE_100 = 20*log10(max(x)*sqrt(256)/(norm(x - x_LMMSE)));
    end
    if i== 4
        cov_X = 1000 * (L.') * L ;    
        x_LMMSE = cov_X*(H.')* (inv(H*cov_X*(H.') + var_W*var_W*identity_)) * y +  zeros(256,1);
        plot(x_LMMSE, 'Color','c','LineWidth',2,'DisplayName','\lambda = 1000');
        hold on;
        result_MMSE_1000 = 20*log10(max(x)*sqrt(256)/(norm(x - x_LMMSE)));
    end
end

title('Linear MMSE Estimated Signals via changing \lambda and Ground Truth Signal');
xlabel('samples');
ylabel('magnitude');
legend();

% PSNR calculations for each estimators

PSNR_dB_ML = 20*log10(max(x)*sqrt(256)/(norm(x - x_ML)));


