/* function declarations */
void InsertAtom(double r[3], int type);

/* functions */
void InsertAtom(double r[3], int type)
{
    atoms[atomNumber].id = atoms[atomNumber - 1].id + 1;
    atoms[atomNumber].r[0] = r[0];
    atoms[atomNumber].r[1] = r[1];
    atoms[atomNumber].r[2] = r[2];
    atoms[atomNumber].type = type;
    atomNumber++;
}
