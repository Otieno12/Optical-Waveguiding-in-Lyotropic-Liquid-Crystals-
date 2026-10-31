// LC_zigzag_alpha_single_frame.ijm  (ImageJ / Fiji macro)
//
// Estimates the zigzag angle alpha (ray path vs. tube axis) from ONE image frame.
// If the open image is a stack, only the currently displayed slice is used.
//
// Steps
//   1. Draw a line along the tube axis -> image is rotated so the axis is horizontal.
//   2. Draw a rectangle around a straight section of tube -> image is cropped.
//   3. AUTOMATIC: intensity-weighted centre across the tube in each column gives a trace y_c(x);
//      a sliding-window line fit gives the local alpha.
//   4. MANUAL (optional): click the zigzag vertices (wall bounces) in order with the Multi-point
//      tool; alpha is computed for every segment between consecutive points.
//
// Usage: open the frame (PNG/TIFF) or scroll a stack to the frame you want, then run the macro.
// Pixels must be square (check Image > Properties) or the angles will be wrong.

// ---------- parameters you may want to change ----------
colStep     = 4;    // columns (px) averaged per sample along the tube
win         = 12;   // samples per local fit (window = win*colStep px); keep below half a zigzag period
blurSigma   = 2;    // Gaussian blur (px) before the automatic analysis
frac        = 0.5;  // threshold fraction between column min and max used for the centroid
minContrast = 8;    // skip columns whose (max - min) intensity is below this
minFill     = 0.8;  // fraction of valid samples needed in a window
// --------------------------------------------------------

if (nImages == 0) exit("Open the image frame first, then run the macro.");

run("Duplicate...", "title=work");          // copies only the current slice
if (bitDepth() != 8) run("8-bit");
workID = getImageID();

// ---- Step 1: tube axis -> rotate so the axis is horizontal ----
setTool("line");
waitForUser("Step 1 of 2",
    "Draw a STRAIGHT LINE along the TUBE AXIS (Straight line tool), then click OK.");
if (selectionType() != 5) exit("No straight line found. Run the macro again.");
getLine(x1, y1, x2, y2, lw);
ang = atan2(y1 - y2, x2 - x1) * 180 / PI;   // counter-clockwise angle of the axis (image y points down)
if (ang > 90)   ang = ang - 180;
if (ang <= -90) ang = ang + 180;
run("Select None");
setBackgroundColor(0, 0, 0);                // black fill for the corners created by the rotation
run("Rotate... ", "angle=" + ang + " grid=1 interpolation=Bicubic enlarge");
print("Tube axis was at " + d2s(ang, 2) + " deg; image rotated so the axis is horizontal.");

// ---- Step 2: region of interest -> crop ----
setTool("rectangle");
waitForUser("Step 2 of 2",
    "Draw a RECTANGLE around a straight section of tube\n" +
    "(full tube width, a little background each side, avoid the black corners), then click OK.");
if (selectionType() != 0) exit("No rectangle found. Run the macro again.");
run("Crop");
cropID = getImageID();
W = getWidth();
H = getHeight();
nC = floor(W / colStep);
stepW = maxOf(1, floor(win / 2));           // window advance (samples)

// ---- Automatic analysis (on a blurred copy) ----
run("Duplicate...", "title=blurred");
run("Gaussian Blur...", "sigma=" + blurSigma);
xs = newArray(nC);
ys = newArray(nC);
for (c = 0; c < nC; c++) {
    xs[c] = c * colStep + colStep / 2;
    ys[c] = NaN;
    prof = newArray(H);
    for (yy = 0; yy < H; yy++) {
        s = 0;
        for (i = 0; i < colStep; i++) s += getPixel(c * colStep + i, yy);
        prof[yy] = s / colStep;
    }
    Array.getStatistics(prof, mn, mx);
    if (mx - mn >= minContrast) {
        thr = mn + frac * (mx - mn);
        sw = 0;
        sy = 0;
        for (yy = 0; yy < H; yy++) {
            if (prof[yy] > thr) {
                wgt = prof[yy] - thr;
                sw += wgt;
                sy += wgt * yy;
            }
        }
        if (sw > 0) ys[c] = sy / sw;
    }
}
close();                                    // blurred copy

Plot.create("Centroid trace", "x (px)", "y (px, image coordinates)", xs, ys);
Plot.show();

Table.create("Zigzag_alpha");
row = 0;
allAlpha = newArray(0);
for (c = 0; c + win <= nC; c += stepW) {
    n = 0; sx = 0; sY = 0; sxy = 0; sxx = 0;
    for (i = c; i < c + win; i++) {
        if (!isNaN(ys[i])) {
            n++;
            sx += xs[i];
            sY += ys[i];
            sxy += xs[i] * ys[i];
            sxx += xs[i] * xs[i];
        }
    }
    if (n >= win * minFill) {
        den = n * sxx - sx * sx;
        slope = -(n * sxy - sx * sY) / den;      // minus: y-up convention
        a = atan(abs(slope)) * 180 / PI;
        Table.set("x_center_px", row, sx / n, "Zigzag_alpha");
        Table.set("slope_yup", row, slope, "Zigzag_alpha");
        Table.set("alpha_deg", row, a, "Zigzag_alpha");
        row++;
        allAlpha = Array.concat(allAlpha, a);
    }
}
Table.update("Zigzag_alpha");

print("---- AUTOMATIC alpha (centroid trace) ----");
if (allAlpha.length == 0) {
    print("No valid windows found. Lower minContrast or win, or check the ROI.");
} else {
    Array.getStatistics(allAlpha, amin, amax, amean, asd);
    sorted = Array.copy(allAlpha);
    Array.sort(sorted);
    amed = sorted[floor(sorted.length / 2)];
    print("windows: " + allAlpha.length);
    print("mean   : " + d2s(amean, 2) + " deg");
    print("median : " + d2s(amed, 2) + " deg");
    print("std    : " + d2s(asd, 2) + " deg");
    print("min/max: " + d2s(amin, 2) + " / " + d2s(amax, 2) + " deg");
}

// ---- Manual mode (optional): click the zigzag vertices ----
selectImage(cropID);
run("Select None");
setTool("multipoint");
waitForUser("Manual mode (optional)",
    "Zoom in and use the Multi-point tool to click the zigzag VERTICES\n" +
    "(the wall bounces) IN ORDER along the tube, at the centre of the bright path.\n" +
    "Need at least 2 points. Click OK when done, or just click OK to skip.");
if (selectionType() == 10) {
    getSelectionCoordinates(px, py);
    if (px.length >= 2) {
        Table.create("Manual_alpha");
        manual = newArray(0);
        for (k = 0; k < px.length - 1; k++) {
            dx = px[k + 1] - px[k];
            dy = py[k + 1] - py[k];
            am = atan2(abs(dy), abs(dx)) * 180 / PI;
            Table.set("segment", k, k + 1, "Manual_alpha");
            Table.set("dx_px", k, dx, "Manual_alpha");
            Table.set("dy_px", k, dy, "Manual_alpha");
            Table.set("alpha_deg", k, am, "Manual_alpha");
            manual = Array.concat(manual, am);
        }
        Table.update("Manual_alpha");
        Array.getStatistics(manual, mmin, mmax, mmean, msd);
        print("---- MANUAL alpha (clicked vertices) ----");
        print("segments: " + manual.length);
        print("mean    : " + d2s(mmean, 2) + " deg");
        print("std     : " + d2s(msd, 2) + " deg");
        print("min/max : " + d2s(mmin, 2) + " / " + d2s(mmax, 2) + " deg");
    }
}

dir = getDirectory("Choose a folder to save the CSV files (Cancel to skip)");
if (dir != "") {
    if (isOpen("Zigzag_alpha")) Table.save(dir + "zigzag_alpha_auto.csv", "Zigzag_alpha");
    if (isOpen("Manual_alpha")) Table.save(dir + "zigzag_alpha_manual.csv", "Manual_alpha");
}
