// ImageJ / Fiji macro: attenuation along a dye-doped LC tube from a waveguiding video
// All lengths in um, using your scale bar (default: 10 um = 75 px -> 0.1333 um/px)
//
// Steps:
//   1. Open the AVI and average all frames
//   2. You draw a straight line along the tube axis, starting at the excitation spot
//   3. You draw a straight line across the tube, wall to wall (diameter)
//   4. Macro extracts the averaged profile, fits I = a*exp(-k*x) + c,
//      and reports 1/e length, dB/um, and length in tube diameters
//
// Run: Plugins > Macros > Run...  (or drag onto Fiji)

// ---------- settings ----------
Dialog.create("Waveguide attenuation (um)");
Dialog.addNumber("Scale bar length (um):", 10);
Dialog.addNumber("Scale bar length (pixels):", 75);
Dialog.addNumber("Averaging strip width across tube (um):", 2.3);
Dialog.addNumber("Fit starts at distance from spot (um):", 5.3);
Dialog.addMessage("Next: choose the AVI file.");
Dialog.show();
barUm = Dialog.getNumber();
barPx = Dialog.getNumber();
stripUm = Dialog.getNumber();
fitFromUm = Dialog.getNumber();

umPerPx = barUm / barPx;                 // 10/75 = 0.1333 um per pixel
stripPx = round(stripUm / umPerPx);      // line width in pixels (about 17)
if (stripPx < 1) stripPx = 1;
unit = "um";

// ---------- open + average ----------
path = File.openDialog("Select the waveguide AVI");
run("AVI...", "open=[" + path + "] use convert");   // use = virtual stack, convert = gray
origTitle = getTitle();
run("Z Project...", "projection=[Average Intensity]");
avgTitle = getTitle();
run("32-bit");
close(origTitle);
selectImage(avgTitle);
run("Set Scale...", "distance=" + barPx + " known=" + barUm + " unit=um");
run("Enhance Contrast", "saturated=0.35");   // display only, data unchanged

// ---------- axis line ----------
setTool("line");
waitForUser("Axis line",
    "Draw a STRAIGHT line along the tube axis:\n" +
    "start at the bright excitation spot, end far along the tube\n" +
    "(about 47-53 um, i.e. 350-400 px). Then click OK.");
if (selectionType() != 5) exit("Please draw a straight line selection.");
getLine(x1, y1, x2, y2, lw);
run("Line Width...", "line=" + stripPx);
profile = getProfile();                        // averaged across the strip
n = profile.length;
distPx = newArray(n);
dist = newArray(n);                            // distance from spot in um
stepPx = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2)) / (n - 1);
for (i = 0; i < n; i++) {
    distPx[i] = i * stepPx;
    dist[i] = distPx[i] * umPerPx;
}
// keep axis drawn on image for the record
setColor(255);
run("Add Selection...");

// ---------- diameter line ----------
run("Line Width...", "line=1");
waitForUser("Diameter line",
    "Draw a straight line ACROSS the tube, from one wall to the other,\n" +
    "perpendicular to the axis. Then click OK.");
if (selectionType() != 5) exit("Please draw a straight line selection.");
getLine(dx1, dy1, dx2, dy2, dlw);
D = sqrt(pow(dx2 - dx1, 2) + pow(dy2 - dy1, 2)) * umPerPx;   // um

// ---------- fit (from fitFromUm onward) ----------
first = 0;
while (first < n && dist[first] < fitFromUm) first++;
m = n - first;
if (m < 15) exit("Too few points after the fit start. Draw a longer axis line or lower the fit start.");
fx = Array.slice(dist, first, n);
fy = Array.slice(profile, first, n);
Fit.doFit("Exponential with Offset", fx, fy);   // y = a*exp(-bx) + c, x in um
a = Fit.p(0);
k = Fit.p(1);                     // 1/um
c = Fit.p(2);
r2 = Fit.rSquared;
Lze = 1 / k;                      // 1/e length, um
dBper = 4.343 * k;                // dB per um
dBperD = dBper * D;               // dB per tube diameter
Lratio = Lze / D;

// ---------- results ----------
print("\\Clear");
print("=== Waveguide attenuation ===");
print("Image: " + path);
print("Scale: " + d2s(umPerPx, 4) + " um/px (" + barUm + " um = " + barPx + " px)");
print("Strip width: " + d2s(stripUm, 2) + " um (" + stripPx + " px)   fit from: " + d2s(fitFromUm, 2) + " um");
print("Fit: I = a*exp(-k*x) + c");
print("  a = " + d2s(a, 2) + "   k = " + d2s(k, 4) + " 1/um   c = " + d2s(c, 2));
print("  R^2 = " + d2s(r2, 3));
print("1/e length        : " + d2s(Lze, 1) + " um");
print("Attenuation       : " + d2s(dBper, 3) + " dB/um   (" + d2s(dBper * 10000, 0) + " dB/cm)");
print("Tube diameter D   : " + d2s(D, 1) + " um");
print("1/e length / D    : " + d2s(Lratio, 2));
print("Loss per diameter : " + d2s(dBperD, 2) + " dB");

// profile + fit plot
fitCurve = newArray(m);
for (i = 0; i < m; i++) fitCurve[i] = a * exp(-k * fx[i]) + c;
Plot.create("Axial profile", "distance from spot (um)", "intensity (a.u.)");
Plot.setColor("black");
Plot.add("circles", dist, profile);
Plot.setColor("red");
Plot.add("line", fx, fitCurve);
Plot.setLegend("data\tsingle exp fit", "top-right");
Plot.show();

// table for export
Table.create("Profile");
Table.setColumn("distance_px", distPx);
Table.setColumn("distance_um", dist);
Table.setColumn("intensity", profile);
Table.update;
// To save: File > Save As > Results... from the "Profile" table window.
