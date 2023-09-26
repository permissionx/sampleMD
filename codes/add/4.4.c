/* classes */
struct Atom
{
    ...
    double minDirection[3];
};

/* function declarations */
void MinDirection_SD();

/* functions declarations */
void MinDirection_SD()
{
    int i, d;
    for (i = 0; i < atomNumber; i++)
    {
        for (d = 0; d < 3; d++)
        {
            atoms[i].minDirection[d] = atoms[i].force[d];
        }
    }
}
