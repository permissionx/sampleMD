/* classes */
struct Atom
{
    ...
    double startR_lineMin[3];
};

/* global variables */
double delta_lineMinimize;

/* function declarations */
void LineMinimize();

/* functions */
void LineMinimize()
{
    double startPotentialEnergy;
    double k, c;
    double slop;
    double lambda_t;
    int i, d;

    k = 0.5;
    c = 0.5;
    slop = 0;
    lambda_t = 100;

    startPotentialEnergy = totalPotentialEnergy;
    for (i = 0; i < atomNumber; i++)
    {
        for (d = 0; d < 3; d++)
        {
            slop += atoms[i].minDirection[d] * atoms[i].minDirection[d];
            atoms[i].startR_lineMin[d] = atoms[i].r[d];
        }
    }
    slop = sqrt(slop);

    lambda_t /= k;
    do
    {
        lambda_t *= k;
        for (i = 0; i < atomNumber; i++)
        {
            for (d = 0; d < 3; d++)
            {
                atoms[i].r[d] = atoms[i].startR_lineMin[d] + lambda_t * atoms[i].minDirection[d];
            }
        }
        PBC_r();
        NeighborList(1);
        Potential(1, 0);
    } while (totalPotentialEnergy > startPotentialEnergy - c * slop * lambda_t);
}
