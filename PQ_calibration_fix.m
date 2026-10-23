% =========================================================================
% PQ_calibration_fix.m
% Fixing the short-calibration weakness of the online integer detector
% (D-TMA + local rule, ideas: Mohsen Karimi) on the low-cost meter data.
% Why 15 s "counting" calibration is weak: at 0.1 false alarms/s only
% ~1.5 alarms happen in 15 s - a rate cannot be estimated from 1-2 events.
% Fix: during calibration the device logs the score (peak/noise) of ALL
% noise peaks (hundreds), fits the upper TAIL of their distribution and
% extrapolates the threshold k for the requested false-alarm rate.
% Methods compared (same 15 s, nothing injected):
%   count  : smallest k with count/T <= budget            (old method)
%   exp    : exponential tail fit of the top peak scores
%   pareto : power-law (Pareto) tail fit  - more cautious for heavy tails
%   par/2  : pareto aimed at half the budget (safety margin)
% Robustness: each 60 s recording -> 4 rotations (calibrate on one 15 s
% block, test on the other ~45 s) -> 16 calibration/test pairs.
% Text output only, no toolbox. Uses the 4 CSV files in MATLAB Drive.
% =========================================================================
clear; clc; tic;
fprintf('==== PQ_calibration_fix  (%s) ====\n', datestr(now));
files = {'V_I_Waveforms_330W_Resistive_Load.csv',          'resistive 330 W';
         'V_I_Waveforms_0.5HP_Induction_Motor_No_Load.csv', 'induction motor';
         'V_I_Waveforms_50W_Ceiling_Fan_Aged.csv',          'ceiling fan 50 W';
         'V_I_Waveforms_SMPS_Mobile_Charger.csv',           'SMPS charger'};
Budg = [0.3 0.1];  nRot = 4;  QS = 2048;  Gc = 4;  kLow = 3;
Amps = [0.3 0.2 0.15];
Wsec = 0.5e-3;  tau1 = 100e-6;  tau2 = 400e-6;
MN = {'count', 'exp', 'pareto', 'par/2'};
nF = size(files,1); nB = numel(Budg); nM = numel(MN); nA = numel(Amps);
FA = nan(nF, nRot, nM, nB);  KK = nan(nF, nRot, nM, nB);
HIT = zeros(nM, nB, nA);  TOT = zeros(nM, nB, nA);

for f = 1:nF
    fn = files{f,1};
    if ~isfile(fn), fprintf('%s missing - run PQ_four_loads.m first\n', fn); continue; end
    [x, fs, f0, bad] = loadUniform(fn);
    N = numel(x);  P = params(fs, Wsec, f0);
    K1 = ceil(P.w1/2);  K2 = ceil(P.w2/2);  m = round(log2(P.L));
    D = (K1-1) + (K2-1);  G = max(P.w1, P.Tmin);  warm = round(0.2*fs);
    [y, s] = streamFilter(round(x*QS), K1, K2, m, Gc, warm);
    e = abs(y);
    % all noise peaks above a low level (what the device logs)
    dC = decide(y, s, kLow, G, warm);
    pC = dC - D;  okC = pC >= 1;  dC = dC(okC);  pC = pC(okC);
    okC = ~bad(pC);  dC = dC(okC);  pC = pC(okC);
    zC = e(dC) ./ s(dC);
    % rotations
    seg = floor((N - warm)/nRot);
    rot = zeros(nRot, 2);
    for r = 1:nRot, rot(r,:) = [warm + (r-1)*seg + 1, warm + r*seg]; end
    % transients injected once per amplitude (whole file, gaps avoided)
    [sh, pk] = shape(fs, tau1, tau2);  rng(500 + f);
    INJ = cell(nA, 2);
    for a = 1:nA
        [xr, ev] = injectOnly(x, bad, fs, Amps(a), sh, pk);
        [yr, sr] = streamFilter(round(xr*QS), K1, K2, m, Gc, warm);
        INJ{a,1} = {yr, sr};  INJ{a,2} = ev;
    end
    fprintf('\n%s: %d logged noise peaks in %.0f s (%.0f /s)\n', files{f,2}, ...
        numel(zC), (N - warm)/fs, numel(zC)/((N - warm)/fs));
    for r = 1:nRot
        a0 = rot(r,1);  b0 = rot(r,2);
        inCal = pC >= a0 & pC <= b0;
        Tc = sum(~bad(a0:b0))/fs;
        tmask = true(N,1); tmask(1:warm) = false; tmask(a0:b0) = false;
        tmask(bad) = false;  Tt = sum(tmask)/fs;
        for bI = 1:nB
            k = calK(zC(inCal), Tc, Budg(bI), kLow);
            KK(f, r, :, bI) = k;
            for mI = 1:nM
                d = decide(y, s, k(mI), G, warm) - D;  d = d(d >= 1);
                FA(f, r, mI, bI) = sum(tmask(d))/Tt;
                for a = 1:nA
                    yr = INJ{a,1}{1};  sr = INJ{a,1}{2};  ev = INJ{a,2};
                    ev = ev(tmask(ev));
                    dr = decide(yr, sr, k(mI), G, warm) - D;  dr = dr(dr >= 1);
                    dr = dr(tmask(dr));
                    h = hitIdx(dr, ev, P.tol);
                    HIT(mI, bI, a) = HIT(mI, bI, a) + sum(h);
                    TOT(mI, bI, a) = TOT(mI, bI, a) + numel(ev);
                end
            end
        end
        fprintf('  rot %d (cal %4.1f-%4.1f s): FA/s on unseen test  [budget 0.3 | 0.1]  ', ...
            r, a0/fs, b0/fs);
        for mI = 1:nM
            fprintf('%s %.2f|%.2f  ', MN{mI}, FA(f,r,mI,1), FA(f,r,mI,2));
        end
        fprintf('\n');
    end
end

fprintf('\n#####################################################################\n');
fprintf('SUMMARY over %d calibration/test pairs (4 loads x %d rotations)\n', nF*nRot, nRot);
for bI = 1:nB
    B = Budg(bI);
    fprintf('\nBudget %.2f false alarms/s\n', B);
    fprintf('  %-7s  FA/s mean  median   max | pairs > budget  > 2x budget | Se ', 'method');
    fprintf('%.2f  ', Amps); fprintf('p.u.\n');
    for mI = 1:nM
        v = FA(:,:,mI,bI); v = v(~isnan(v));
        fprintf('  %-7s   %6.3f  %6.3f  %6.3f |      %2d/%-2d        %2d/%-2d    |', ...
            MN{mI}, mean(v), myq(v,0.5), max(v), sum(v > B), numel(v), ...
            sum(v > 2*B), numel(v));
        for a = 1:nA, fprintf(' %.3f', HIT(mI,bI,a)/max(1,TOT(mI,bI,a))); end
        fprintf('\n');
    end
end
fprintf('\nNote: with a correct calibration, ~half of the pairs can still exceed\n');
fprintf('the budget by chance (Poisson). Look at mean and "> 2x budget".\n');
fprintf('\nTotal run time: %.1f s\n==== END ====\n', toc);

% =========================================================================
%                              FUNCTIONS
% =========================================================================
function k = calK(z, T, B, kLow)
% thresholds from the logged peak scores z (T seconds, nothing injected)
z  = sort(z(:), 'descend');  n = numel(z);  k = kLow*ones(1,4);
if n == 0, return; end
% count
nb = floor(B*T);
if nb + 1 <= n, k(1) = z(nb + 1); else, k(1) = kLow; end
% tail: top 10 % of peaks, at least 30, at most 300
n0 = min(n, max(30, min(300, round(0.1*n))));
u0 = z(n0);  R0 = n0/T;  top = z(1:n0);
targets = [B B B/2];
for j = 1:3
    Bt = targets(j);
    if R0 <= Bt, k(j+1) = max(u0, kLow); continue; end
    if j == 1
        lam = 1/max(mean(top - u0), eps);
        k(2) = u0 + log(R0/Bt)/lam;
    else
        al = 1/max(mean(log(top/u0)), eps);
        k(j+1) = u0*(R0/Bt)^(1/al);
    end
end
k = max(k, kLow);
end

function [y, s, mx] = streamFilter(q, K1, K2, m, Gc, warm)
% genuine sample-by-sample causal loop, integer values only
N = numel(q); y = zeros(N,1); s = ones(N,1);
b1 = zeros(K1,1); c1 = zeros(K1,1); b2 = zeros(K2,1); c2 = zeros(K2,1);
Dl = K1 - 1; dl = zeros(max(Dl,1),1);
S1 = 0; S2 = 0; T1 = 0; T2 = 0; i1 = 1; i2 = 1; id = 1;
acc = 0; sumA = 0; mx = 0; M2 = 2^m;
for n = 1:N
    xn = q(n);
    S1 = S1 + xn - b1(i1);  b1(i1) = xn;          % running sum 1
    S2 = S2 + S1 - c1(i1);  c1(i1) = S1;          % running sum 2 (TMA*K1^2)
    i1 = i1 + 1; if i1 > K1, i1 = 1; end
    if Dl > 0, xd = dl(id); dl(id) = xn; id = id + 1; if id > Dl, id = 1; end
    else, xd = xn; end
    h  = K1*K1*xd - S2;                           % exact DC rejection
    T1 = T1 + h - b2(i2);  b2(i2) = h;
    T2 = T2 + T1 - c2(i2); c2(i2) = T1;
    i2 = i2 + 1; if i2 > K2, i2 = 1; end
    e = abs(T2);
    if abs(S2) > mx, mx = abs(S2); end
    if abs(h)  > mx, mx = abs(h);  end
    if e       > mx, mx = e;       end
    if n <= warm                                  % warm-up: plain mean |y|
        sumA = sumA + e;  acc = floor(sumA*M2/n);
    else
        sc  = floor(acc/M2);                      % = acc >> m
        acc = acc + min(e, Gc*sc) - sc;           % winsorised EMA
    end
    if acc > mx, mx = acc; end
    y(n) = T2;  s(n) = max(floor(acc/M2), 1);
end
end

function [det, tem] = decide(y, s, k, G, warm)
% same result as the per-sample loop: local max above k*s, chain merge
% within G samples (keep larger e/s), emit when G samples pass without a
% stronger peak. tem = emission sample (decision time).
e = abs(y); n = numel(e); i = (2:n-1)';
c = i(e(i) > e(i-1) & e(i) >= e(i+1) & e(i) > k*s(i) & i > warm);
det = zeros(numel(c),1); tem = det; nd = 0; pend = 0; ps = 0;
for j = 1:numel(c)
    i = c(j); sc = e(i)/s(i);
    if pend > 0 && i - pend < G
        if sc > ps, pend = i; ps = sc; end
    else
        if pend > 0, nd = nd + 1; det(nd) = pend; tem(nd) = pend + G; end
        pend = i; ps = sc;
    end
end
if pend > 0, nd = nd + 1; det(nd) = pend; tem(nd) = pend + G; end
det = det(1:nd); tem = tem(1:nd);
end

function [h, jdx] = hitIdx(det, ev, tol)
det = det(:); used = false(size(det)); h = false(numel(ev),1); jdx = ones(numel(ev),1);
for j = 1:numel(ev)
    if isempty(det), break; end
    d = abs(det - ev(j)); d(used) = inf; [dm, k] = min(d);
    if dm <= tol, used(k) = true; h(j) = true; jdx(j) = k; end
end
end

function [x, ev] = injectOnly(bg, bseg, fs, A0, sh, pk)
N = numel(bg); ev = []; p = round(0.02*fs); mg = round(0.01*fs);
while true
    p = p + round((0.005 + 0.115*rand)*fs);
    if p > N - round(0.02*fs) - numel(sh), break; end
    if ~any(bseg(max(1,p-mg) : min(N, p+mg+numel(sh)))), ev(end+1, 1) = p; end %#ok<AGROW>
end
s = zeros(N, 1);
for j = 1:numel(ev)
    a  = (2*(rand > 0.5) - 1) * (0.6 + 0.8*rand) * A0;
    ii = (ev(j)-pk+1):(ev(j)-pk+numel(sh));
    s(ii) = s(ii) + a*sh;
end
x = bg(:) + s;
end

function [x, fs, f0, bad, ng] = loadUniform(fn)
T = readtable(fn, 'VariableNamingRule', 'preserve');
t = double(T.time_s(:)); v = double(T.voltage_V(:));
ok = ~isnan(t) & ~isnan(v); t = t(ok); v = v(ok);
[t, ia] = unique(t); v = v(ia);
dt = diff(t); md = median(dt); fs = round(1/md);
tu = (t(1) : 1/fs : t(end))';
x  = interp1(t, v, tu, 'linear');
g  = find(dt > 2*md); ng = numel(g);
bad = false(numel(tu), 1); pad = 0.02;
for k = 1:ng
    i1 = max(1, floor((t(g(k)) - pad - t(1))*fs) + 1);
    i2 = min(numel(tu), ceil((t(g(k)+1) + pad - t(1))*fs) + 1);
    bad(i1:i2) = true;
end
x = x - mean(x);
[tc, nUp] = upcross(x); f0 = (nUp-1)/(tc(end)-tc(1))*fs;
x = x/(sqrt(2)*sqrt(mean(x.^2)));
end

function P = params(fs, Wsec, f0)
P.fs = fs; P.f0 = f0;
P.W  = max(3, round(Wsec*fs));
P.w1 = max(3, round(3.3*P.W));
P.w2 = max(3, round(0.06*P.W));
P.L  = round(10/f0*fs);          % 10 fundamental cycles (IEC 61000-4-30 window)
P.Tmin = round(2e-3*fs);         % physical minimum interval between transients
P.tol  = max(2, round(0.5e-3*fs));
P.N0 = round(fs/f0);             % one cycle (for cycle-difference comparator)
P.k  = sqrt(2*log(P.L));        % threshold factor (default, not calibrated)
P.kw = NaN;                      % NaN -> universal threshold for db4
end

function q = myq(x, p)          % own quantile (no toolbox)
x = sort(x(:)); n = numel(x);
if n == 0, q = NaN; return; end
pos = 1 + (n-1)*p; lo = floor(pos); hi = ceil(pos);
q = x(lo) + (pos-lo)*(x(hi) - x(lo));
end

function [sh, pk] = shape(fs, tau1, tau2)
tt = (0:round(10*tau2*fs))'/fs;
sh = exp(-tt/tau2) - exp(-tt/tau1);
sh = sh / max(sh);
[~, pk] = max(sh);
end

function [tc, nUp] = upcross(x)  % up-crossings with hysteresis, linear interpolation
a = 0.1*max(abs(x)); st = 0; tc = zeros(ceil(numel(x)/2), 1); nUp = 0; last = 2;
for i = 2:numel(x)
    if x(i-1) <= 0 && x(i) > 0, last = i; end
    if st <= 0 && x(i) > a
        if st < 0
            k = last;
            nUp = nUp + 1;
            tc(nUp) = (k-1) + (-x(k-1))/(x(k) - x(k-1));
        end
        st = 1;
    elseif st >= 0 && x(i) < -a
        st = -1;
    end
end
tc = tc(1:nUp);
end
