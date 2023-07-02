/* constants */
#define MAX_NEIGHBOR_NUMBER 2000 // maximum number of neighbors

/* classes */
struct AtomNeighbor
{
    int index;
    double distance;
    double dr[3];  // from atom to neighbor
};

struct Atom
{
    int id, type;
    double r[3];
    double reR[3];
    double boxReR[3];
    int neighborNumber;
    struct AtomNeighbor neighbors[MAX_NEIGHBOR_NUMBER];
};

/* global variables */
double neighborCutoff;

/* function declarations */
void ConstructNeighborList();

/* functions */
void ConstructNeighborList()
{
    int i, j;
    int d;
    double distance;
    double dr[3];
    int neighborIndex_i, neighborIndex_j;

    for (i = 0; i < atomNumber; i++)
    {
        atoms[i].neighborNumber = 0;
    }
    for (i = 0; i < atomNumber; i++)
    {
        for (j = i + 1; j < atomNumber; j++)
        {
            PBC_dr(i, j, dr);
            distance = sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]);
            if (distance < neighborCutoff)
            {
                neighborIndex_i = atoms[i].neighborNumber;
                neighborIndex_j = atoms[j].neighborNumber;
                atoms[i].neighbors[neighborIndex_i].index = j;
                atoms[j].neighbors[neighborIndex_j].index = i;
                atoms[i].neighbors[neighborIndex_i].distance = distance;
                atoms[j].neighbors[neighborIndex_j].distance = distance;
                for (d = 0; d < 3; d++)
                {
                    atoms[i].neighbors[neighborIndex_i].dr[d] = dr[d];
                    atoms[j].neighbors[neighborIndex_j].dr[d] = -dr[d];
                }
                atoms[i].neighborNumber++;
                atoms[j].neighborNumber++;
            }
        }
    }
}

/* main */
int main()
{
    /* parameters */
    neighborCutoff = 1.1;

    /* processing*/
    ConstructStdCrystal_BCC(1.0, 5);
    ConstructNeighborList();

    /* output */
    int n, d, index;
    printf("count index neighbor_id x y z dx dy dz dr\n");
    for (n = 0; n < atoms[0].neighborNumber; n++)
    {
        index = atoms[0].neighbors[n].index;
        printf("%d %d %d ", n+1, index, atoms[index].id);
        for (d = 0; d < 3; d++)
        {
            printf("%f ", atoms[index].r[d]);
        }
        for (d = 0; d < 3; d++)
        {
            printf("%f ", atoms[0].neighbors[n].dr[d]);
        }
        printf("%f \n", atoms[0].neighbors[n].distance);
    }

    return 0;
}
