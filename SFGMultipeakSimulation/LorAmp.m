function Amp = LorAmp(A,Omega,gamma,IR)
% 仅适用于 JCP,124,114705-1,2006 的模拟
% 注意不同文献中公式的+-会有区别
    Amp = A ./ (Omega - IR - 1i .* gamma);
end
