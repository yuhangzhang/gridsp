
#include "../nifticlib-2.0.0/niftilib/nifti1_io.h"

#include "../vxl-1.14.0/core/vgui/vgui.h"
#include "../vxl-1.14.0/core/vgui/vgui_viewer2D_tableau.h"
#include "../vxl-1.14.0/core/vgui/vgui_image_tableau.h"

#include "../vxl-1.14.0/core/vil/vil_image_view.h"
#include "../vxl-1.14.0/core/vil/vil_load.h"
#include "../vxl-1.14.0/core/vil/vil_convert.h"

#include "../supervoxel/superVoxel.h"

#include "gridSV.h"

//exe height_of_volumne size_of_superpixels

int main(int argc, char * argv[])
{
	nifti_image * nim =  nifti_image_read("C:\\cygwin64\\home\\Yuhang\\visceral\\visceral-dataset\\trainingset\\volumes\\temp\\10000014_1_CT_wb.nii",1);

	printf("ndim=%d\n",nim->ndim);

	for(int i=0;i<=nim->ndim;i++)
	{
		printf("%d\n",nim->dim[i]);
	}

	printf("nvox=%d\n",nim->nvox);
	printf("nbyper=%d\n",nim->nbyper);
	printf("datatype=%d\n",nim->datatype);
	printf("cal_min=%f\n",nim->cal_min);
	printf("cal_max=%f\n",nim->cal_max);
	//getchar();

	float * mydata = (float*)nim->data;

	float minvox,maxvox;

	minvox = mydata[0];
	maxvox = mydata[0];

	for(int i=0; i<nim->nvox ;i++)
	{
		if(minvox>=mydata[i]) minvox=mydata[i];
		if(maxvox<=mydata[i]) maxvox=mydata[i];
	}

	float maxgap = maxvox-minvox;


	vil_image_view<float> im(nim->dim[1],nim->dim[2],atoi(argv[1]));

	for(int i=0;i<im.ni();i++)
	{
		for(int j=0;j<im.nj();j++)
		{
			for(int k=0;k<im.nplanes();k++)
			{
				im(i,j,k) = (mydata[k*im.ni()*im.nj()+i*im.nj()+j+510*im.ni()*im.nj()]-minvox)/maxgap;
				//if(i%20==1) im(i,j,k)=im(i-1,j,k);
				//else im(i,j,k)=1;
			}
		}
	}

	//vil_image_view<vxl_byte> tim = vil_load("other.jpg");
	//
	//vil_convert_cast(tim,im);

	//printf("here %d %d\n",im.ni(),im.nj());

	//for(int i=0;i<im.ni();i++)
	//{
	//	for(int j=0;j<im.nj();j++)
	//	{
	//		im(i,j) = im(i,j)/255.0;
	//	}
	//}

	//getchar();

	//superVoxel<float> sv(im,atoi(argv[2]));
	gridSV<float> sv(im,atoi(argv[2]));

	printf("here2\n");

	vil_image_view<vxl_byte> im2(im.ni(),im.nj(),3);

	for(int i=0;i<im2.ni()-1;i++)
	{
		for(int j=0;j<im2.nj()-1;j++)
		{
			for(int k=0;k<3;k++)
			{
				im2(i,j,k) = im(i,j,0)*255;				
			}
			if(sv.getLabel(i,j,0)!=sv.getLabel(i+1,j,0)
				||sv.getLabel(i,j,0)!=sv.getLabel(i,j+1,0))
			{
				im2(i,j,0)=200;
				im2(i,j,1)=0;
				im2(i,j,2)=0;
			}
			//printf("%d ",sv.getLabel(i,j,0));

		}
		//printf("\n");
	}


	vgui::init(argc,argv);

	vgui_image_tableau_new image(im2);
    vgui_viewer2D_tableau_new viewer(image);
    vgui::run(viewer, image->width(), image->height());
	


	nifti_image_free( nim );


	return 0;
}

