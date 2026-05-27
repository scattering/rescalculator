int Res3_Init() ;
void Res3_shutdown() ;

void Res3_setDblprec( int prec ) ;
char *Res3_getErrmsg() ;
int Res3_GetInput(char *name, char **data) ;
int Res3_GetOutput(char *name, char **data) ;
int Res3_GetNextInput(char **name, char **data, int init) ;
int Res3_GetNextOutput(char **name, char **data, int init) ;
int Res3_SetOp(int argc, char **argv) ;
void Res3_SetDopt(double h, double k, double l, double e) ;
Conv3_resinfo *Res3_NewResinfo() ;
Conv3_conv3 *Res3_NewConv3s() ;
void Res3_CopyResinfo(Conv3_resinfo *srcRes, Conv3_resinfo *destRes) ;
int Res3_CopyConvDatapt(Conv3_conv3 *s, int is, Conv3_conv3 *d, int id) ;
void Res3_FreeConv3s(Conv3_conv3 *lPtr) ;
int Res3_GetConvNpts(Conv3_conv3 *lPtr) ;
int Res3_checkConvNpts(int np, Conv3_conv3 *conv3pt) ;
int Res3_Conv3Prep(Conv3_conv3 *conv3pt) ;
Conv3_conv3 *Res3_Conv3MeshPrep(Conv3_conv3 *conv3pt,
				Conv3_conv3 *interconv, int mesh) ;
void Res3_unloadResinfo(Conv3_resinfo *respt) ;
void Res3_loadResinfo(Conv3_resinfo *respt) ;
