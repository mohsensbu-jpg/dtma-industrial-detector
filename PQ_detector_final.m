% =========================================================================
% PQ_detector_final.m
% FINAL online integer event detector (core ideas: Mohsen Karimi):
%   1) D-TMA band-pass, integer running sums, one parameter W
%   2) local noise scale = streaming 99 % quantile of |y|
%      (integer; up 99 steps / down 1 step, step = q >> 12)
%   3) alarm = local peak > k * q99, peaks merged within Tmin
%   4) 15 s calibration: log all noise peaks, fit the Pareto tail of
%      their scores, choose k for the requested false-alarm rate
% This script, on the low-cost meter (resistive load):
%   - calibrates on 15 s, measures false alarms on the unseen rest
%   - injects 0.3 p.u. transients (~2/s) and measures detection + delay
%   - prints cost, memory and the largest internal integer
%   - saves the 3 portfolio figures (PNG) in MATLAB Drive
% Needs V_I_Waveforms_330W_Resistive_Load.csv (run PQ_four_loads.m once)
% Text output only, no toolbox.
% =========================================================================
clear; clc; tic;
fprintf('==== PQ_detector_final  (%s) ====\n', datestr(now));
fn = 'V_I_Waveforms_330W_Resistive_Load.csv';
if ~isfile(fn), error('file missing - run PQ_four_loads.m first'); end
QS = 2048;  Tcal = 15;  Budget = 0.1;  kLow = 1;  A0 = 0.3;
Wsec = 0.5e-3;  tau1 = 100e-6;  tau2 = 400e-6;
cBlue = [0.12 0.38 0.70];  cRed = [0.80 0.20 0.15];
cGrey = [0.55 0.55 0.55];  cGreen = [0.10 0.55 0.30];

%% ---------------- detector set-up and calibration -----------------------
[x, fs, f0, bad] = loadUniform(fn);
N = numel(x);  P = params(fs, Wsec, f0);
K1 = ceil(P.w1/2);  K2 = ceil(P.w2/2);
D = (K1-1) + (K2-1);  G = max(P.w1, P.Tmin);  warm = round(0.2*fs);
fprintf('fs=%d Hz  W=%.2f ms  K1=%d K2=%d  delay D=%.2f ms  Tmin=%.1f ms\n', ...
    fs, Wsec*1e3, K1, K2, D/fs*1e3, G/fs*1e3);
[y, q, mx0] = streamQ(round(x*QS), K1, K2, warm);
e = abs(y);
calm = false(N,1);  calm(warm+1 : warm+round(Tcal*fs)) = true;
tmask = ~calm & ~bad;  tmask(1:warm) = false;
dC = decide(y, q, kLow, G, warm);  pC = dC - D;
ok = pC >= 1;  dC = dC(ok);  pC = pC(ok);
ok = calm(pC) & ~bad(pC);
z = e(dC(ok)) ./ q(dC(ok));
k = parK(z, sum(calm & ~bad)/fs, Budget, kLow);
fprintf('calibration: %d noise peaks logged in %d s -> k = %.2f (budget %.2f /s)\n', ...
    numel(z), Tcal, k, Budget);
d = decide(y, q, k, G, warm) - D;  d = d(d >= 1);  d = d(tmask(d));
Tt = sum(tmask)/fs;
fprintf('false alarms on unseen %.1f s (nothing injected): %d = %.3f /s\n', Tt, numel(d), numel(d)/Tt);

%% ---------------- injected transients -----------------------------------
[sh, pk] = shape(fs, tau1, tau2);  rng(7);
[xr, ev] = injectRate(x, bad | calm, fs, A0, sh, pk, 0.2, 0.8);
ev = ev(tmask(ev));
[yr, qr, mx1] = streamQ(round(xr*QS), K1, K2, warm);
[dO, tO] = decide(yr, qr, k, G, warm);
pos = dO - D;  keep = pos >= 1;  pos = pos(keep);  tO = tO(keep);
keep = tmask(pos);  pos = pos(keep);  tO = tO(keep);
[h, jdx] = hitIdx(pos, ev, P.tol);
lat = (tO(jdx(h)) - ev(h))/fs*1e3;
fprintf('injected %.1f p.u. transients: %d, detected %d (%.1f %%), extra alarms %d\n', ...
    A0, numel(ev), sum(h), 100*mean(h), numel(pos) - sum(h));
fprintf('decision delay: median %.1f ms, max %.1f ms\n', myq(lat, 0.5), max(lat));
mx = max(mx0, mx1);
fprintf('largest |internal integer|: %d (int32 limit 2147483647 -> %s)\n', mx, ...
    ternary(mx < 2^31-1, 'OK', 'OVERFLOW'));
fprintf('per sample: filter 9 add/sub + 1 multiply; quantile tracker 1 compare,\n');
fprintf('  1 shift, 1 add; decision 3 compares + 1 multiply (k*q)\n');
fprintf('memory: %d integers of buffers + 8 state variables\n', 2*K1 + 2*K2 + K1 - 1);

%% ---------------- figure 1: detection example ---------------------------
Lw = round(0.6*fs);  st = 0;
for j = 1:numel(ev)
    a = ev(j) - round(0.05*fs);  b = a + Lw - 1;
    if a < warm + 1 || b > N - D, continue; end
    if any(bad(a:b)) || any(calm(a:b)), continue; end
    if sum(ev >= a & ev <= b) >= 3, st = a; break; end
end
if st == 0, st = find(tmask, 1); end
w  = st:(st + Lw - 1);
tm = (w - st)/fs*1e3;
evw = ev(ev >= w(1) & ev <= w(end));
dw  = pos(pos >= w(1) & pos <= w(end));
f1 = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 900 620]);
subplot(3,1,1);
plot(tm, xr(w), 'Color', cBlue, 'LineWidth', 0.8); hold on;
plot((evw - st)/fs*1e3, xr(evw), 'v', 'Color', cRed, 'MarkerFaceColor', cRed, 'MarkerSize', 6);
ylabel('voltage (p.u.)'); grid on; box off;
title('1) Real grid voltage from a low-cost meter + transients (red marks)');
subplot(3,1,2);
plot(tm, abs(yr(w + D))/QS, 'Color', cBlue, 'LineWidth', 0.8); hold on;
plot(tm, k*qr(w + D)/QS, '--', 'Color', cRed, 'LineWidth', 1.2);
ylabel('|filter output|'); grid on; box off;
legend({'D-TMA output (integer, online)', 'adaptive threshold'}, ...
    'Location', 'northeast', 'Box', 'off');
title('2) The filter removes the 50 Hz wave; the threshold follows the local noise');
subplot(3,1,3);
plot(tm, zeros(size(tm)), 'Color', cGrey); hold on;
h1 = plot((evw - st)/fs*1e3, ones(size(evw)), 'o', 'Color', cRed, 'MarkerSize', 9, 'LineWidth', 1.5);
h2 = plot((dw - st)/fs*1e3, ones(size(dw)), 'x', 'Color', cGreen, 'MarkerSize', 10, 'LineWidth', 2);
ylim([-0.5 1.8]); yticks([]); xlabel('time (ms)'); box off;
legend([h1 h2], {'true transient', 'alarm'}, 'Location', 'northeast', ...
    'Box', 'off', 'Orientation', 'horizontal');
title(sprintf('3) Alarms, %.0f ms after each event', myq(lat, 0.5)));
print(f1, 'fig1_detection_example.png', '-dpng', '-r130');  close(f1);

%% ---------------- figure 2: detection vs amplitude ----------------------
% PQ_adaptive_threshold.m, method static/q99: 4 loads x 4 calibration
% windows, 11.1 min of unseen test data, ~2 transients/s
A    = [0.2 0.3 0.5];
Se1  = [0.592 0.910 0.997];     % budget 0.10 /s, measured 0.077 /s
Se2  = [0.352 0.778 0.988];     % budget 0.03 /s, measured 0.014 /s
f2 = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 760 480]);
hold on;
p1 = plot(A, 100*Se1, '-o', 'Color', cBlue, 'LineWidth', 2.5, 'MarkerSize', 8, 'MarkerFaceColor', cBlue);
p2 = plot(A, 100*Se2, '-s', 'Color', cGrey, 'LineWidth', 2, 'MarkerSize', 7, 'MarkerFaceColor', cGrey);
for i = 1:3
    text(A(i), 100*Se1(i) + 5, sprintf('%.0f%%', 100*Se1(i)), ...
        'HorizontalAlignment', 'center', 'Color', cBlue, 'FontWeight', 'bold');
end
xlim([0.17 0.53]); ylim([0 108]); grid on; box off;
xticks(A); xlabel('transient amplitude (p.u. of nominal peak)');
ylabel('transients detected (%)');
legend([p1 p2], {'false alarms 0.08 /s (measured)', 'false alarms 0.014 /s (measured)'}, ...
    'Location', 'southeast', 'Box', 'off');
title({'Low-cost Raspberry-Pi meter, 10 kHz, real noise, 4 loads', ...
       'online integer detector, 15 s calibration, 11 min of unseen data'});
print(f2, 'fig2_sensitivity.png', '-dpng', '-r130');  close(f2);

%% ---------------- figure 3: block diagram -------------------------------
f3 = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1000 360]);
ax = axes('Position', [0 0 1 1]); axis(ax, [0 100 0 36]); axis off; hold on;
bx = {'Sensor', {'isolated PT / CT', 'accelerometer', '4-20 mA, 0-10 V'};
      'ADC', {'timer-triggered', 'DMA, 10-50 kHz', 'no timing jitter'};
      'D-TMA filter', {'integer running sums', '~10 ops / sample', 'one parameter W'};
      'Local noise', {'streaming 99 %', 'quantile of |y|', 'shift update'};
      'Decision', {'peak > k \cdot noise', 'merge within T_{min}', 'delay ~3 ms'};
      'Outputs', {'relay / LED', 'Modbus RTU / TCP', 'event log'}};
nb = size(bx, 1);  bw = 13.5;  gap = (100 - nb*bw)/(nb + 1);
for i = 1:nb
    x0 = gap + (i-1)*(bw + gap);
    if i >= 3 && i <= 5, fc = [0.90 0.94 1.00]; else, fc = [0.96 0.96 0.96]; end
    rectangle('Position', [x0 12 bw 18], 'Curvature', 0.12, ...
        'FaceColor', fc, 'EdgeColor', cBlue, 'LineWidth', 1.3);
    text(x0 + bw/2, 26.5, bx{i,1}, 'HorizontalAlignment', 'center', ...
        'FontWeight', 'bold', 'FontSize', 11);
    for j = 1:3
        text(x0 + bw/2, 22 - (j-1)*3.4, bx{i,2}{j}, ...
            'HorizontalAlignment', 'center', 'FontSize', 8.5);
    end
    if i < nb
        annotation('arrow', [(x0 + bw + 0.4)/100, (x0 + bw + gap - 0.4)/100], ...
            [21/36 21/36], 'Color', cBlue, 'LineWidth', 1.3);
    end
end
text(50, 6.5, ['Set-up: W (event width), T_{min} (minimum interval), ' ...
    'false-alarm budget  +  15 s automatic calibration (noise-peak tail fit)'], ...
    'HorizontalAlignment', 'center', 'FontSize', 10, 'Color', [0.2 0.2 0.2]);
text(50, 33.5, 'Low-cost online event detector (shaded = runs on the microcontroller)', ...
    'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
print(f3, 'fig3_block_diagram.png', '-dpng', '-r130');  close(f3);

fl = {'fig1_detection_example.png', 'fig2_sensitivity.png', 'fig3_block_diagram.png'};
for i = 1:3
    dd = dir(fl{i});
    if isempty(dd), fprintf('NOT saved: %s\n', fl{i});
    else, fprintf('saved: %-28s %6.1f kB\n', fl{i}, dd.bytes/1024); end
end
fprintf('fig 1 window: %d transients, %d alarms\n', numel(evw), numel(dw));
fprintf('Total run time: %.1f s\n==== END ====\n', toc);

% =========================================================================
%                              FUNCTIONS
% =========================================================================
function [y, q, mx] = streamQ(qin, K1, K2, warm)
% THE DETECTOR CORE: integer, causal, sample by sample
N = numel(qin); y = zeros(N,1); q = ones(N,1);
b1 = zeros(K1,1); c1 = zeros(K1,1); b2 = zeros(K2,1); c2 = zeros(K2,1);
Dl = K1 - 1; dl = zeros(max(Dl,1),1);
S1 = 0; S2 = 0; T1 = 0; T2 = 0; i1 = 1; i2 = 1; id = 1;
sumA = 0; qq = 0; mx = 0;
for n = 1:N
    xn = qin(n);
    S1 = S1 + xn - b1(i1);  b1(i1) = xn;          % running sum 1
    S2 = S2 + S1 - c1(i1);  c1(i1) = S1;          % running sum 2
    i1 = i1 + 1; if i1 > K1, i1 = 1; end
    if Dl > 0, xd = dl(id); dl(id) = xn; id = id + 1; if id > Dl, id = 1; end
    else, xd = xn; end
    h  = K1*K1*xd - S2;                           % exact DC rejection
    T1 = T1 + h - b2(i2);  b2(i2) = h;
    T2 = T2 + T1 - c2(i2); c2(i2) = T1;
    i2 = i2 + 1; if i2 > K2, i2 = 1; end
    e = abs(T2);
    if n <= warm                                  % start: 3 x mean |y|
        sumA = sumA + e;  qq = 3*floor(sumA/n);
    else                                          % streaming 99 % quantile
        stp = max(floor(qq/4096), 1);             % = qq >> 12
        if e > qq, qq = qq + 99*stp; else, qq = qq - stp; end
    end
    mx = max([mx, abs(S2), abs(h), e, qq]);
    y(n) = T2;  q(n) = max(qq, 1);
end
end

function k = parK(z, T, B, kLow)
% Pareto (power-law) tail fit of logged noise-peak scores -> threshold
z = sort(z(:), 'descend');  n = numel(z);
if n < 30 || T <= 0, k = 7*kLow; return; end
n0 = min(n, max(30, min(300, round(0.1*n))));
u0 = z(n0);  R0 = n0/T;  top = z(1:n0);
if R0 <= B, k = max(u0, kLow); return; end
al = 1/max(mean(log(top/u0)), eps);
k  = max(u0*(R0/B)^(1/al), kLow);
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

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
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
