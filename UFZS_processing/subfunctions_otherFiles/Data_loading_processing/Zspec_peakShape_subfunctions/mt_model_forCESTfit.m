function out1=mt_model_forCESTfit(t2a,ra,rb,mb0,r,t2b,sigtau,delta_MT,delta,w1,lineshape)

% function out1=mtmodel(t2a,ra,rb,mb0,r,t2b,delta,w1,lineshape)
% Note: may need to fix  t2a , ra, and rb to get good fit;
% Scott Swanson August 13 2024 (c)

% DK notes: 
%   Most of these equations are found in Henkelman et al, MRM 1993
%   delta and delta_MT should be in Hz
%   w1 should be in rad/s
%   sigtau [unitless] is for Kubo-Tomita lineshape
%   Updated 11/11/25 to include sigtau for Kubo-Tomita lineshape
%   Updated 5/14/26 for lineshape options: {'lorentzian','gaussian','superlorentzian','kubo-tomita'}


% Fixing values for white matter fitting
% t2a=0.065;
% ra=0.4;
% rb=1;

if nargin<10 %shift over variables if sigtau not included
    w1=delta;
    delta=delta_MT;
    delta_MT=sigtau;
    sigtau=1;
end

if ~exist('lineshape','var')
    warning('Lineshape function not specified! Setting to Gaussian by default...')
    lineshape='gaussian';
end

rrfa=t2a.*w1.^2./(1+(t2a.*2*pi*delta).^2);  %eq10, Henkelman paper
                                            %NOTE: typo in paper! Missing
                                            %square of parentheses term

% Switch calculation of parameters based upon variable lineshape
switch lineshape
    case 'lorentzian'   % Lorentzian semisolid lineshape

        rrfb=t2b.*w1.^2./(1+(t2b.*2*pi*(delta-delta_MT)).^2);

    case 'gaussian'     % Gaussian semisolid lineshape -- DK note: makes curve "lumpier"

        rrfb=t2b*sqrt(pi/2)*w1.^2.*exp(-((2*pi*(delta-delta_MT)*t2b).^2)/2);

    case 'superlorentzian' % Super-Lorentzian semisolid lineshape

        rrfb=RF_superlorentzian(t2b,delta_MT,w1,delta);

    case 'kubo-tomita'  % Kubo-Tomita semisolid lineshape

        rrfb=RF_KuboTomita(t2b,sigtau,delta_MT,w1,delta);

end
                                            

%These three below come from eq9, Henkelman paper by factoring the
%numerator and multiplying out the denominator
numerator = r.*(ra+mb0*rb)+ra.*(rb+rrfb); 
denom1 = (ra+rrfa).*(rb+rrfb);
denom2 = r.*(ra+rrfa+mb0.*(rb+rrfb));   

mtspec=numerator./(denom1 + denom2);

% plot(numerator,'linewidth',2)
% hold on
% % plot(denom1,'linewidth',2)
% % plot(denom2,'linewidth',2)
% % plot(denom3,'linewidth',2)
% % legend('num','den1','den2','den3')
% hold off

%plot(delta,mtspec,'o')
out1 = mtspec;
