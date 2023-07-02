/* classes */ 
struct Atom
{
    ...
    double lastForce_CG[3];
};

/* function declarations */
void MinDirection_CG();

/* functions */
void MinDirection_CG()
{
    int i, d;
    double beta_up, beta_down, beta;
    double directionCheck;
    if (nStep == 0)
    {
        for (i = 0; i < atomNumber; i++)
        {
            for (d = 0; d < 3; d++)
            {
                atoms[i].minDirection[d] = atoms[i].force[d];
                atoms[i].lastForce_CG[d] = atoms[i].force[d];
            }
        }
    }
    else
    {
        beta_up = 0;
        beta_down = 0;
        for (i = 0; i < atomNumber; i++)
        {
            beta_up += VecDotMul(atoms[i].force, atoms[i].force);
            beta_down += VecDotMul(atoms[i].lastForce_CG, atoms[i].lastForce_CG);
        }
        beta = beta_up / beta_down;
        directionCheck = 0;
        for (i = 0; i < atomNumber; i++)
        {
            for (d = 0; d < 3; d++)
            {
                atoms[i].minDirection[d] = atoms[i].force[d] + beta * atoms[i].minDirection[d];
                atoms[i].lastForce_CG[d] = atoms[i].force[d];
                directionCheck += atoms[i].minDirection[d] * atoms[i].force[d];
            }
        }
        if (directionCheck < 0)
        {
            for (i = 0; i < atomNumber; i++)
            {
                for (d = 0; d < 3; d++)
                {
                    atoms[i].minDirection[d] = atoms[i].force[d];
                }
            }
        }
    }
}
