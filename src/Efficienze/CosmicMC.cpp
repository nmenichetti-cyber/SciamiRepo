#include <iostream>
#include <cmath>
#include <iostream>
#include <cmath>


//runna con MC(1000000, 3.5)
class ray {
    double cosTheta;
    double phi;
    double x, y;
    double exp;
    double size_x;
    double size_y;

public:
    ray(double exp, double size_x, double size_y) 
        : exp(exp), size_x(size_x), size_y(size_y) {
        Throw();
    }

    void Throw() {
        phi = gRandom->Uniform() * 2 * M_PI;
        cosTheta = pow(gRandom->Uniform(), 1.0 / (exp + 1)); 
        x = gRandom->Uniform(size_x); 
        y = gRandom->Uniform(size_y);
    }

    double CosTheta() const { return cosTheta; }
    double Phi() const { return phi; }
    double X() const { return x; }
    double Y() const { return y; }

    double X(double z) const {
        return x - sqrt(1 - pow(CosTheta(), 2)) / CosTheta() * cos(Phi()) * z;
    }

    double Y(double z) const {
        return y - sqrt(1 - pow(CosTheta(), 2)) / CosTheta() * sin(Phi()) * z;
    }
};

class detector {
    double xmin, xmax, ymin, ymax, z;
public:
    detector(double xmin, double xmax, double ymin, double ymax, double z)
        : xmin(xmin), xmax(xmax), ymin(ymin), ymax(ymax), z(z) {}

    bool Check(const ray& cr) {
        double x_on = cr.X(z);
        double y_on = cr.Y(z);
        return x_on > xmin && x_on < xmax && y_on > ymin && y_on < ymax;
    }
};

void MC(int n, double exponent) {
    ray r(exponent, .50, .28);
    detector down(0., .50, 0., .40, -0.205);
    detector up(0., .50, 0., .40, 0.205);

    TH2D* hUp = new TH2D("hUp",
        "Intercette detector UP; x [m]; y [m]",
        100, 0., .50,
        100, 0., .28);

    TH2D* hDown = new TH2D("hDown",
        "Intercette detector DOWN; x [m]; y [m]",
        100, 0., .50,
        100, 0., .28);

    int pairs = 0, triples = 0;

    for(int i = 0; i < n; i++) {
        r.Throw(); 
        pairs   += up.Check(r);
        triples += down.Check(r) && up.Check(r);

            if (up.Check(r)) {
            hUp->Fill(r.X(+.27), r.Y(+.27));
        }

        if (down.Check(r)) {
            hDown->Fill(r.X(-.37), r.Y(-.37));
        }
    }

    double A = (double)triples / pairs;

    TCanvas* c = new TCanvas("c", "Intercette", 1000, 400);
    c->Divide(2,1);

    c->cd(1);
    hUp->Draw("COLZ");

    c->cd(2);
    hDown->Draw("COLZ");

    cout << "Pairs = " << pairs << endl;
    cout << "Triples =" << triples << endl;
    cout << "A = " << A <<endl;
}