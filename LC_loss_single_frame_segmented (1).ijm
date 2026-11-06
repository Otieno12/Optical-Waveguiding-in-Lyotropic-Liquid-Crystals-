// LC_loss_single_frame_segmented.ijm  (ImageJ / Fiji macro)
//
// Estimates the attenuation along a (possibly bent or wavy) tube from ONE image frame.
// You trace the tube axis with a SEGMENTED line; the intensity is averaged across the tube width
// along that path, and an exponential is fitted:
//
//      I(z) - bg = I0 * exp(-gamma * z)
//
//   gamma            attenuation coefficient (1/um)
//   L = 1/gamma      1/e propagation length (um)
//   loss (dB/cm)     = 4.343 * gamma(1/cm)
//
// Steps
//   1. Segmented line tool: click along the tube axis (from just after the entry spot to the far
//      end), double-click to finish.
//   2. Rectangle on a dark background region (away from the tube and its halo).
//   3. Enter the scale (default: 120 px = 10 um), tube diameter (optional), line width, smoothing.
//   4. Look at the ln(I - bg) plot, then enter the distance range for the fit.
//
// Use the RAW, linear frame (no gamma / contrast changes), without saturated pixels if possible.
// If the image is a stack, the currently displayed slice is used. Needs a reasonably recent Fiji.

// ---------- parameters you may want to change ----------
minNetFrac = 0.03;   // ignore samples whose (I - bg) is below this fraction of the maximum
// --------------------------------------------------------

if (nImages == 0) exit("Open the image frame first, then run the macro.");

run("Duplicate...", "title=lossimg");        // only the current slice
if (bitDepth() == 24) run("8-bit");
run("Set Scale...", "distance=0 known=0 unit=pixel");   // work in pixels; the dialog converts to um
imgID = getImageID();
satLevel = 255;
if (bitDepth() == 16) satLevel = 65535;

// ---- Step 1: tube axis (segmented line) ----
setTool("polyline");
waitForUser("Step 1 of 2",
    "Use the SEGMENTED LINE tool: click along the TUBE AXIS from just after the entry spot\n" +
    "to the far end, DOUBLE-CLICK to finish, then click OK.");
st = selectionType();
if (st != 5 && st != 6 && st != 7) exit("No line selection found. Run the macro again.");
getSelectionCoordinates(xp, yp);

// ---- Step 2: background ----
setTool("rectangle");
waitForUser("Step 2 of 2",
    "Draw a RECTANGLE on a DARK BACKGROUND region (away from the tube and its glow),\n" +
    "then click OK.");
if (selectionType() != 0) exit("No rectangle found. Run the macro again.");
getStatistics(areaBg, bg);
print("Background level: " + d2s(bg, 2));

// ---- Parameters ----
Dialog.create("Loss measurement");
Dialog.addNumber("Scale: length in pixels:", 120);
Dialog.addNumber("Scale: equals (um):", 10);
Dialog.addNumber("Tube diameter (um; 0 = unknown):", 0);
Dialog.addNumber("Line width for averaging across the tube (px):", 40);
Dialog.addNumber("Smoothing along the tube (px):", 5);
Dialog.show();
scalePx = Dialog.getNumber();
scaleUm = Dialog.getNumber();
tubeD   = Dialog.getNumber();
lw      = Dialog.getNumber();
smoothPx = Dialog.getNumber();
pixSize = scaleUm / scalePx;                 // um per pixel
print("Pixel size: " + d2s(pixSize, 5) + " um/px");

// ---- Profile along the segmented line ----
selectImage(imgID);
makeSelection("polyline", xp, yp);
Roi.setStrokeWidth(lw);
lenPx = getValue("Length");                  // path length in px (image is uncalibrated)
prof = getProfile();
N = prof.length;
dpx = lenPx / (N - 1);                       // px between samples (about 1)
sm = round(smoothPx / dpx);
sp = newArray(N);
dist = newArray(N);
for (i = 0; i < N; i++) {
    s = 0;
    cnt = 0;
    for (j = i - sm; j <= i + sm; j++) {
        if (j >= 0 && j < N) {
            s += prof[j];
            cnt++;
        }
    }
    sp[i] = s / cnt;
    dist[i] = i * dpx * pixSize;
}
nSat = 0;
for (i = 0; i < N; i++) if (prof[i] >= satLevel) nSat++;
if (nSat > 0) print("WARNING: " + nSat + " profile samples are saturated; start the fit after them.");
print("Path length: " + d2s(dist[N - 1], 1) + " um");

Array.getStatistics(sp, spMin, spMax);
minNet = minNetFrac * (spMax - bg);
net = newArray(N);
lnI = newArray(N);
peakIdx = 0;
lastIdx = 0;
for (i = 0; i < N; i++) {
    net[i] = sp[i] - bg;
    lnI[i] = NaN;
    if (net[i] > minNet) {
        lnI[i] = log(net[i]);                // natural log
        lastIdx = i;
    }
    if (sp[i] > sp[peakIdx]) peakIdx = i;
}

Plot.create("Profile: ln(I - bg) vs distance", "distance (um)", "ln(I - bg)", dist, lnI);
Plot.show();

// ---- Fit range ----
Dialog.create("Fit range");
Dialog.addMessage("Look at the plot and choose the clean, straight, decaying part.\n" +
    "Start after any saturated or entry region. Distances in um.");
Dialog.addNumber("Fit start:", dist[peakIdx]);
Dialog.addNumber("Fit end:", dist[lastIdx]);
Dialog.show();
s0 = Dialog.getNumber();
s1 = Dialog.getNumber();

// ---- Fit ----
n = 0; sx = 0; sy = 0; sxy = 0; sxx = 0; syy = 0;
for (i = 0; i < N; i++) {
    if (dist[i] >= s0 && dist[i] <= s1 && !isNaN(lnI[i])) {
        n++;
        sx += dist[i];
        sy += lnI[i];
        sxy += dist[i] * lnI[i];
        sxx += dist[i] * dist[i];
        syy += lnI[i] * lnI[i];
    }
}
if (n < 5) exit("Fewer than 5 valid samples in the fit range. Widen the range.");
den = n * sxx - sx * sx;
slope = (n * sxy - sx * sy) / den;
icpt = (sy - slope * sx) / n;
r = (n * sxy - sx * sy) / sqrt(den * (n * syy - sy * sy));
if (slope >= 0) exit("The fitted slope is not negative (no decay in this range).");
gamma = -slope;                              // 1/um
Lprop = 1 / gamma;                           // um

// sensitivity of gamma to the background level
gLow  = fitGamma(dist, sp, 0.5 * bg, s0, s1, minNetFrac);
gHigh = fitGamma(dist, sp, 1.5 * bg, s0, s1, minNetFrac);

// ---- Plot data + fit ----
xf = newArray(2);
yf = newArray(2);
xf[0] = s0; xf[1] = s1;
yf[0] = icpt + slope * s0;
yf[1] = icpt + slope * s1;
Plot.create("Fit: ln(I - bg) vs distance", "distance (um)", "ln(I - bg)", dist, lnI);
Plot.setColor("red");
Plot.setLineWidth(2);
Plot.add("line", xf, yf);
Plot.show();

// ---- Report ----
print("---- Attenuation from a single frame ----");
print("Fit range      : " + d2s(s0, 2) + " to " + d2s(s1, 2) + " um (" + n + " samples)");
print("R^2            : " + d2s(r * r, 3));
print("gamma          : " + d2s(gamma, 4) + " per um");
print("L (1/e)        : " + d2s(Lprop, 2) + " um");
print("loss           : " + d2s(4.343 * gamma * 10000, 0) + " dB/cm");
if (tubeD > 0) print("L / D          : " + d2s(Lprop / tubeD, 2));
if (!isNaN(gLow) && !isNaN(gHigh)) {
    print("Background sensitivity (bg x0.5 / x1.5):");
    print("   L (1/e)     : " + d2s(1 / gLow, 2) + " / " + d2s(1 / gHigh, 2) + " um");
    print("   loss        : " + d2s(4.343 * gLow * 10000, 0) + " / " + d2s(4.343 * gHigh * 10000, 0) + " dB/cm");
}
print("Note: the decay of the dye glow tracks the guided excitation light and includes dye absorption.");

// ---- Save profile ----
Table.create("Loss_profile");
for (i = 0; i < N; i++) {
    Table.set("distance_um", i, dist[i], "Loss_profile");
    Table.set("intensity_smoothed", i, sp[i], "Loss_profile");
    Table.set("net_intensity", i, net[i], "Loss_profile");
    Table.set("ln_net", i, lnI[i], "Loss_profile");
}
Table.update("Loss_profile");
dir = getDirectory("Choose a folder to save loss_profile.csv (Cancel to skip)");
if (dir != "") Table.save(dir + "loss_profile.csv", "Loss_profile");


// ---------------- helper ----------------
function fitGamma(d, intensity, bgLevel, a, b, frac) {
    nn = d.length;
    mx = 0;
    for (k = 0; k < nn; k++) if (intensity[k] - bgLevel > mx) mx = intensity[k] - bgLevel;
    thr = frac * mx;
    m = 0; qx = 0; qy = 0; qxy = 0; qxx = 0;
    for (k = 0; k < nn; k++) {
        v = intensity[k] - bgLevel;
        if (d[k] >= a && d[k] <= b && v > thr) {
            lv = log(v);
            m++;
            qx += d[k];
            qy += lv;
            qxy += d[k] * lv;
            qxx += d[k] * d[k];
        }
    }
    if (m < 5) return NaN;
    dd = m * qxx - qx * qx;
    sl = (m * qxy - qx * qy) / dd;
    if (sl >= 0) return NaN;
    return -sl;
}
