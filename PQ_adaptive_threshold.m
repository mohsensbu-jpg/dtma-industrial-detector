% =========================================================================
% PQ_adaptive_threshold.m
% Online integer D-TMA detector (ideas: Mohsen Karimi) on the low-cost
% meter: making the false-alarm rate hold on unseen data.
% Finding so far: a one-time 15 s calibration fails when the meter noise
% changes during operation. Fix = extend his "local" idea from the noise
% SCALE to the THRESHOLD itself:
%   static  : Pareto tail fit on the first 15 s, then fixed
%   rolling : every 1 s, re-fit the tail on the noise peaks of the last
%             15 s (peaks that raised an alarm are excluded), causal
% Two local noise scales:
%   mean    : recursive mean |y| (as before)
%   q99     : streaming 99 % quantile of |y| (integer, shift update) -
%             follows the tail of the noise directly
%   rollRare: rolling k learned from noise only (real transients are rare,
%             a few per hour, so they do not move the threshold)
% Realistic event rate: ~2 transients/s, none during the first 15 s.
% 4 loads x 4 rotations = 16 calibration/test pairs (~12 min of test).
% Text output only, no toolbox. Uses the 4 CSV files in MATLAB Drive.
% =========================================================================
clear; clc; tic;
fprintf('==== PQ_adaptive_threshold  (%s) ====\n', datestr(now));
files = {'V_I_Waveforms_330W_Resistive_Load.csv',          'resistive 330 W';
         'V_I_Waveforms_0.5HP_Induction_Motor_No_Load.csv', 'induction motor';
         'V_I_Waveforms_50W_Ceiling_Fan_Aged.csv',          'ceiling fan 50 W';
         'V_I_Waveforms_SMPS_Mobile_Charger.csv',           'SMPS charger'};
Budg = [0.1 0.03];  nRot = 4;  QS = 2048;  Gc = 4;  kLowS = [3 1];  Twin = 15;
Amps = [0.5 0.3 0.2];
Wsec = 0.5e-3;  tau1 = 100e-6;  tau2 = 400e-6;
MN = {'static/mean', 'rolling/mean', 'rollRare/mean', 'static/q99', 'rolling/q99', 'rollRare/q99'};
nF = size(files,1); nB = numel(Budg); nM = 6; nA = numel(Amps);
FAn = zeros(nF*nRot, nM, nB);  FAt = zeros(nF*nRot, 1);
HIT = zeros(nM, nB, nA);  TOT = zeros(nM, nB, nA);  pr = 0;

for f = 1:nF
    fn = files{f,1};
    if ~isfile(fn), fprintf('%s missing - run PQ_four_loads.m first\n', fn); continue; end
    [x, fs, f0, bad] = loadUniform(fn);
    N = numel(x);  P = params(fs, Wsec, f0);
    K1 = ceil(P.w1/2);  K2 = ceil(P.w2/2);  m = round(log2(P.L));
    D = (K1-1) + (K2-1);  G = max(P.w1, P.Tmin);  warm = round(0.2*fs);
    [sh, pk] = shape(fs, tau1, tau2);
    [y0, sM0, sQ0] = streamFilter2(round(x*QS), K1, K2, m, Gc, warm);
    seg = floor((N - warm)/nRot);
    fprintf('\n%s\n', files{f,2});
    for r = 1:nRot
        pr = pr + 1;
        a0 = warm + (r-1)*seg + 1;  b0 = warm + r*seg;
        calm = false(N,1);  calm(a0:b0) = true;
        tmask = ~calm & ~bad;  tmask(1:warm) = false;  FAt(pr) = sum(tmask)/fs;
        ordF = [a0:N, warm+1:a0-1]';            % causal order, circular
        posOf = zeros(N,1);  posOf(ordF) = 1:numel(ordF);
        nCal = b0 - a0 + 1;
        cb = [0; cumsum(~bad(ordF))];           % clean samples along ordF
        rng(700 + 10*f + r);
        KVR = cell(2, nB);                      % rolling k from the noise-only run
        for run = 0:nA
            if run == 0
                yy = y0;  SS = {sM0, sQ0};  ev = [];
            else
                [xr, ev] = injectRate(x, bad | calm, fs, Amps(run), sh, pk, 0.2, 0.8);
                [yy, sm, sq] = streamFilter2(round(xr*QS), K1, K2, m, Gc, warm);
                SS = {sm, sq};  ev = ev(tmask(ev));
            end
            e = abs(yy);
            for sc = 1:2
                s = SS{sc};  kLow = kLowS(sc);
                dC = decide(yy, s, kLow, G, warm);
                pC = dC - D;  ok = pC >= 1;  dC = dC(ok);  pC = pC(ok);
                ok = ~bad(pC);  dC = dC(ok);  pC = pC(ok);
                z = e(dC)./s(dC);  zp = posOf(pC);
                Tcal = sum(calm & ~bad)/fs;
                for bI = 1:nB
                    B = Budg(bI);
                    k0 = parK(z(calm(pC)), Tcal, B, kLow);
                    if isnan(k0), k0 = 7*kLow; end
                    % ---- static ----
                    kv = k0*ones(N,1);
                    % ---- rolling ----
                    kvR = k0*ones(N,1);  kAt = inf(size(z));  kprev = k0;
                    st = nCal + 1;
                    while st <= numel(ordF)
                        en = min(numel(ordF), st + fs - 1);
                        w1 = max(1, st - Twin*fs);
                        sel = zp >= w1 & zp < st & z <= kAt;
                        Tw = (cb(st) - cb(w1))/fs;
                        kj = parK(z(sel), Tw, B, kLow);
                        if isnan(kj), kj = kprev; end
                        idx = ordF(st:en);  dIdx = idx + D;  dIdx = dIdx(dIdx <= N);
                        kvR(dIdx) = kj;
                        inb = zp >= st & zp <= en;  kAt(inb) = kj;
                        kprev = kj;  st = en + 1;
                    end
                    if run == 0, KVR{sc, bI} = kvR; end
                    for v = 1:3
                        if v == 1, kk = kv; elseif v == 2, kk = kvR; else, kk = KVR{sc, bI}; end
                        mI = (sc-1)*3 + v;
                        d = decideV(yy, s, kk, G, warm) - D;  d = d(d >= 1);
                        d = d(tmask(d));
                        if run == 0
                            FAn(pr, mI, bI) = numel(d);
                        else
                            h = hitIdx(d, ev, P.tol);
                            HIT(mI, bI, run) = HIT(mI, bI, run) + sum(h);
                            TOT(mI, bI, run) = TOT(mI, bI, run) + numel(ev);
                        end
                    end
                end
            end
        end
        fprintf('  rot %d: FA/s on unseen test, budget %.2f |', r, Budg(1));
        fprintf(' %s %.3f', MN{1}, FAn(pr,1,1)/FAt(pr));
        for mI = 2:nM, fprintf(' | %s %.3f', MN{mI}, FAn(pr,mI,1)/FAt(pr)); end
        fprintf('\n');
    end
end

fprintf('\n#####################################################################\n');
fprintf('SUMMARY: %d calibration/test pairs, %.1f min of unseen test data\n', pr, sum(FAt)/60);
for bI = 1:nB
    B = Budg(bI);
    fprintf('\nBudget %.2f false alarms/s (= %.0f per hour)\n', B, B*3600);
    fprintf('  %-13s pooled FA/s   (per hour)  worst pair  pairs>2x | Se ', 'method');
    fprintf('%.1f  ', Amps); fprintf('p.u.\n');
    for mI = 1:nM
        rate = FAn(:,mI,bI)./FAt;
        pooled = sum(FAn(:,mI,bI))/sum(FAt);
        fprintf('  %-13s    %6.3f      (%5.0f)     %6.3f    %2d/%d  |', MN{mI}, pooled, ...
            pooled*3600, max(rate), sum(rate > 2*B), pr);
        for a = 1:nA, fprintf(' %.3f', HIT(mI,bI,a)/max(1,TOT(mI,bI,a))); end
        fprintf('\n');
    end
end
fprintf('\nTotal run time: %.1f s\n==== END ====\n', toc);

% =========================================================================
%                              FUNCTIONS
% =========================================================================
function k = parK(z, T, B, kLow)
% Pareto (power-law) tail fit of logged noise-peak scores -> threshold
z = sort(z(:), 'descend');  n = numel(z);
if n < 30 || T <= 0, k = NaN; return; end
n0 = min(n, max(30, min(300, round(0.1*n))));
u0 = z(n0);  R0 = n0/T;  top = z(1:n0);
if R0 <= B, k = max(u0, kLow); return; end
al = 1/max(mean(log(top/u0)), eps);
k  = max(u0*(R0/B)^(1/al), kLow);
end

function [y, s, q] = streamFilter2(qin, K1, K2, m, Gc, warm)
% integer causal D-TMA + two local scales: mean |y| (EMA) and the
% streaming 99 % quantile of |y| (up 99 steps / down 1 step, step = q>>12)
N = numel(qin); y = zeros(N,1); s = ones(N,1); q = ones(N,1);
b1 = zeros(K1,1); c1 = zeros(K1,1); b2 = zeros(K2,1); c2 = zeros(K2,1);
Dl = K1 - 1; dl = zeros(max(Dl,1),1);
S1 = 0; S2 = 0; T1 = 0; T2 = 0; i1 = 1; i2 = 1; id = 1;
acc = 0; sumA = 0; M2 = 2^m; qq = 0;
for n = 1:N
    xn = qin(n);
    S1 = S1 + xn - b1(i1);  b1(i1) = xn;
    S2 = S2 + S1 - c1(i1);  c1(i1) = S1;
    i1 = i1 + 1; if i1 > K1, i1 = 1; end
    if Dl > 0, xd = dl(id); dl(id) = xn; id = id + 1; if id > Dl, id = 1; end
    else, xd = xn; end
    h  = K1*K1*xd - S2;
    T1 = T1 + h - b2(i2);  b2(i2) = h;
    T2 = T2 + T1 - c2(i2); c2(i2) = T1;
    i2 = i2 + 1; if i2 > K2, i2 = 1; end
    e = abs(T2);
    if n <= warm
        sumA = sumA + e;  acc = floor(sumA*M2/n);
        qq = 3*floor(acc/M2);
    else
        sc  = floor(acc/M2);
        acc = acc + min(e, Gc*sc) - sc;
        stp = max(floor(qq/4096), 1);
        if e > qq, qq = qq + 99*stp; else, qq = qq - stp; end
    end
    y(n) = T2;  s(n) = max(floor(acc/M2), 1);  q(n) = max(qq, 1);
end
end

function det = decideV(y, s, kv, G, warm)
% as decide(), but with a per-sample threshold factor kv
e = abs(y); n = numel(e); i = (2:n-1)';
c = i(e(i) > e(i-1) & e(i) >= e(i+1) & e(i) > kv(i).*s(i) & i > warm);
det = zeros(numel(c),1); nd = 0; pend = 0; ps = 0;
for j = 1:numel(c)
    i = c(j); sc = e(i)/s(i);
    if pend > 0 && i - pend < G
        if sc > ps, pend = i; ps = sc; end
    else
        if pend > 0, nd = nd + 1; det(nd) = pend; end
        pend = i; ps = sc;
    end
end
if pend > 0, nd = nd + 1; det(nd) = pend; end
det = det(1:nd);
end

function [x, ev] = injectRate(bg, bseg, fs, A0, sh, pk, tmin, tmax)
N = numel(bg); ev = []; p = round(0.05*fs); mg = round(0.01*fs);
while true
    p = p + round((tmin + (tmax - tmin)*rand)*fs);
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
