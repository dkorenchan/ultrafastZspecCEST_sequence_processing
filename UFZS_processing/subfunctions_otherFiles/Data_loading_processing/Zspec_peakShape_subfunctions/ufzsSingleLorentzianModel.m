function Yhat_Single = ufzsSingleLorentzianModel(p,x)
if length(p) < 4, p(4) = 0; end
num = sqrt((p(2)/2)^2 + (x-p(3)).^2) * (p(2)/2);
den = (p(2)/2)^2 + (x-p(3)).^2;
ph  = exp(-1i*p(4));
Yhat_Single = p(1)*real(num./den.*ph);
Yhat_Single = Yhat_Single - min(real(Yhat_Single));
end
