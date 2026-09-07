function Yhat_Single = ufzsSinglePseudoVoigtModel(p,x)
if length(p) < 6, p(6) = 0; end
sigma = p(3) / 2 / sqrt(2*log(2)) * p(4);
G     = exp(-(x-p(5)).^2 ./ 2 ./ sigma^2);
Lnum  = sqrt((p(3)/2)^2 + (x-p(5)).^2) * p(3)/2;
Lden  = (x-p(5)).^2 + (p(3)/2)^2;
Lph   = exp(-1i*(atan((x-p(5))/(p(3)/2)) + p(6)));
L     = Lnum./Lden.*Lph;
Yhat_Single = p(1)*real(p(2)*G + (1-p(2))*L);
end
