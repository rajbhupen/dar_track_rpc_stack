// ROOT Macro to simulate 10-layer RPC detector and visualize muon trajectory
#include <TCanvas.h>
#include <TGraph.h>
#include <TLine.h>
#include <TRandom3.h>
#include <vector>
#include <iostream>

void rpc_detector_sim() {

    const int nLayers = 10;          // Number of RPC layers
    const double layerGap = 5.0;   // Distance between layers in cm
    const int nStrips = 18;         // Number of strips per layer
    const double stripWidth = 3.0;  // Width of each strip in cm
    const double stripGap = 0.5;   // Gap between strips in cm

    // Canvas setup
   
    double scale = 20;
    double totalWidth = nStrips*stripWidth + (nStrips-1)*stripGap;
    double totalHeight = (nLayers-1)*layerGap;
     double xAxisRange = totalWidth * 1.; // Add extra space for labels
    double yAxisRange = totalHeight * 1.; // Add extra space for labels

    TCanvas *c1 = new TCanvas("c1", "Muon Veto Wall Geometry", scale * totalWidth, scale * totalHeight);
    c1->SetGrid();
    c1->Range(-5.0 ,-5.0, xAxisRange + 5.0 , yAxisRange+5.0);




    double inefficiencyRate = 0.1;
    double noiseRate = 0.05;

    // Generate random start and end points for the muon trajectory
    double lxStart = gRandom->Uniform(20, 30); // Detector spans -15 to 15 cm in x
    double lxEnd = gRandom->Uniform(0, 10);
    
    double lyEnd = (nLayers-1) * layerGap +5.0; // Top layer
    double lyStart = 0.0 - 5.0;                    // Bottom layer

    // Draw muon trajectory
    TLine *muonLine = new TLine(lxStart, lyStart, lxEnd, lyEnd);
    muonLine->SetLineColor(kRed);
    muonLine->SetLineWidth(2);

    // Container to store strip intersections
       std::vector<TGraph*> hits;

    // Draw layers and strips
     for (int i = 0; i < nLayers; ++i) {
       cout<<"i: "<<i<<endl;
       double yLayer = i * layerGap;

        for (int j = 0; j < nStrips; ++j) {
	  double sxStart =  j * (stripWidth + stripGap);
            double sxEnd = sxStart + stripWidth;
	    double syStart = yLayer;
	    double syEnd = yLayer;
	    
            // Draw strips
            TLine *stripLine = new TLine(sxStart,syStart,sxEnd,syEnd);
            stripLine->SetLineColor(kBlack);
            stripLine->Draw();
	    stripLine->SetLineWidth(2);

	    
	    //Line
	    double magl = sqrt(pow(lxStart-lxEnd,2.) + pow(lyStart-lyEnd,2.));
	    double dx = (lxEnd - lxStart) / magl;
	    double dy = (lyEnd - lyStart) / magl;
	    
	    //Edge
	    double ex = sxEnd-sxStart;
	    double ey = syEnd-syStart;
	    double mags = sqrt(ex*ex+ey*ey);
	    ex/=mags;
	    ey/=mags;

	    /// Determinant
	    double det = dx * ey - dy * ex;

	    // Check if lines are parallel
	    if (fabs(det) < 1e-9) {
	      std::cout << "Lines are parallel, no intersection." << std::endl;
	      continue; // Skip this edge
	    }

	    // Solve for t and s
	    double t = ((sxStart - lxStart) * ey - (syStart - lyStart) * ex) / det;
	    double s = ((lxStart - sxStart) * dy - (lyStart - syStart) * dx) / det;


	    // Compute intersection point
	    double xIntersect = lxStart + t * dx;
	    double yIntersect = lyStart + t * dy;

	    std::cout << "t = " << t << ", s = " << s << ""<<" xin "<< xIntersect<<" yin "<<yIntersect<< std::endl;

	    
	    if(i==3){
	      cout<<"XXXXXXXXXXXXXXXx"<<endl;
	      cout<<j<<" "<<xIntersect<<" "<<sxStart<<" "<<sxEnd<<endl;
	    }
	    
	    // Check if intersection is valid
	    if (t >= 0 && xIntersect>sxStart && xIntersect<sxEnd) {

	      if (gRandom->Uniform(0, 1) < inefficiencyRate) {
		std::cout << "Strip " << j << " skipped due to inefficiency." << std::endl;
		continue; // Skip this strip
	      }

	      // Draw the intersection point
	      TMarker *intersectionMarker = new TMarker(sxStart+0.5*stripWidth, yIntersect, kFullCircle);
	      intersectionMarker->SetMarkerColor(kRed);
	      intersectionMarker->SetMarkerSize(1.5);
	      intersectionMarker->Draw();
	      
	      std::cout << "Intersection detected at (" << j<<" "<<xIntersect << ", " << yIntersect << ")" << std::endl;


	      //multiplicity
	      // Center of the strip
	      double xStripCenter = sxStart + 0.5 * stripWidth;

	      // Check whether the intersection is on the left or right side
	      bool isRightSide = (xIntersect > xStripCenter);

	      // Generate multiplicity (randomly 1 or 2)
	      int multiplicity = gRandom->Uniform(1, 3); // Random number: 1 or 2 
	      cout<<"multiplicity: "<<multiplicity<<endl;
	      // Place multiplicity markers
	      for (int m = 1; m < multiplicity; ++m) {
		// Determine the neighboring strip based on position
		double xNeighborStripStart = isRightSide
		  ? sxStart + m * (stripWidth + stripGap) // Right neighbor
		  : sxStart - m * (stripWidth + stripGap); // Left neighbor

		// Ensure the neighbor strip index is within valid bounds
		if (xNeighborStripStart < 0 || xNeighborStripStart > nStrips * (stripWidth + stripGap)) {
		  continue; // Skip out-of-bound strips
		}

		// Calculate the marker position for the neighboring strip
		double xNeighborStripCenter = xNeighborStripStart + 0.5 * stripWidth;

		// Draw multiplicity marker
		TMarker *multiplicityMarker = new TMarker(xNeighborStripCenter, yIntersect, kFullCircle);
		multiplicityMarker->SetMarkerColor(kRed); // Blue for multiplicity markers
		multiplicityMarker->SetMarkerSize(1.5);
		multiplicityMarker->Draw();

		std::cout << "Multiplicity marker on strip at x = " << xNeighborStripCenter << ", y = " << yIntersect << std::endl;
	      }

	    }
	    else{
	      

	    }
	    
	      if (gRandom->Uniform(0, 1) < noiseRate) {
		std::cout << "Strip " <<j << " skipped due to inefficiency." << std::endl;
			      // Draw the intersection point
	      TMarker *noiseMarker = new TMarker(sxStart+0.5*stripWidth, syStart, kFullCircle);
	      noiseMarker->SetMarkerColor(kBlack);
	      noiseMarker->SetMarkerSize(1.5);
	      noiseMarker->Draw();
	      
	    }

         }
     }

    muonLine->Draw();

 
double textXOffset = 20.0; // Position text 20 units to the right of the detector's width

for (int i = 0; i < nLayers; ++i) {
    // Calculate the y-position of the current layer
    double yLayer = i * layerGap;

    // Create the layer number text
    TText *layerText = new TText(totalWidth+2., yLayer, Form("%d", i ));
    layerText->SetTextSize(0.03); // Adjust text size
    layerText->SetTextColor(kBlack); // Text color
    layerText->SetTextAlign(22); // Center alignment
    layerText->Draw();
}
    
}
