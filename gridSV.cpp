#include "gridSV.h"

#ifndef _GRIDSV_CPP_
#define _GRIDSV_CPP_

template<class T>
double gridSV<T>::cost(double v1, double v2)
{
	double tempf=-fabs(v1-v2);
	tempf = exp(tempf/0.04)+0.01;
	
	return tempf;
}

template<class T>
gridSV<T>::gridSV(vil_image_view<T> &im, int voxelsize)
{
	_width=im.ni();
	_height=im.nj();
	_layers=im.nplanes();

	volumeElimination vex(voxelsize,voxelsize,voxelsize);
	volumeElimination vey(voxelsize,voxelsize,voxelsize);
	volumeElimination vez(voxelsize,voxelsize,voxelsize);

	volumeElimination::vector3i xlabel;
	volumeElimination::vector3i ylabel;
	volumeElimination::vector3i zlabel;

	xlabel.resize(voxelsize);
	ylabel.resize(voxelsize);
	zlabel.resize(voxelsize);
	_label.resize(_width);

	for(int i=0;i<voxelsize;i++)
	{
		xlabel[i].resize(voxelsize);
		ylabel[i].resize(voxelsize);
		zlabel[i].resize(voxelsize);

		for(int j=0;j<voxelsize;j++)
		{
			xlabel[i][j].resize(voxelsize);
			ylabel[i][j].resize(voxelsize);
			zlabel[i][j].resize(voxelsize);
		}
	}

	for(int i=0;i<_width;i++)
	{
		_label[i].resize(_height);

		for(int j=0;j<_height;j++)
		{
			_label[i][j].resize(_layers);
		}
	}


	for(int i=0,xsign=-1;i<_width;i+=voxelsize)
	{
//printf("i=%d\n",i);
		xsign=-xsign;		
		for(int j=0,ysign=-1;j<_height;j+=voxelsize)
		{
//printf("j=%d\n",j);
			ysign=-ysign;
			for(int k=0,zsign=-1;k<_layers;k+=voxelsize)
			{
//printf("ijk=%d %d %d\n",i,j,k);
				zsign=-zsign;
				for(int ii=0;ii<voxelsize&&i+ii<_width;ii++)
				{
					for(int jj=0;jj<voxelsize&&j+jj<_height;jj++)
					{
						for(int kk=0;kk<voxelsize&&k+kk<_layers;kk++)
						{
							if(ii<voxelsize-1&&i+ii<_width-1)
							{
								double tempd = cost(im(i+ii,j+jj,k+kk),im(i+ii+1,j+jj,k+kk));

								vex.addDataterm(ii,jj,kk,tempd);
								vex.addDataterm(ii+1,jj,kk,tempd);			
								vex.addEdgeterm(ii,jj,kk,ii+1,jj,kk,tempd*-2);

								vey.addDataterm(ii,jj,kk,tempd);
								vey.addDataterm(ii+1,jj,kk,tempd);			
								vey.addEdgeterm(ii,jj,kk,ii+1,jj,kk,tempd*-2);

								vez.addDataterm(ii,jj,kk,tempd);
								vez.addDataterm(ii+1,jj,kk,tempd);			
								vez.addEdgeterm(ii,jj,kk,ii+1,jj,kk,tempd*-2);
							}

							if(jj<voxelsize-1&&j+jj<_height-1)
							{
								double tempd = cost(im(i+ii,j+jj,k+kk),im(i+ii,j+jj+1,k+kk));

								vex.addDataterm(ii,jj,kk,tempd);
								vex.addDataterm(ii,jj+1,kk,tempd);			
								vex.addEdgeterm(ii,jj,kk,ii,jj+1,kk,tempd*-2);

								vey.addDataterm(ii,jj,kk,tempd);
								vey.addDataterm(ii,jj+1,kk,tempd);			
								vey.addEdgeterm(ii,jj,kk,ii,jj+1,kk,tempd*-2);

								vez.addDataterm(ii,jj,kk,tempd);
								vez.addDataterm(ii,jj+1,kk,tempd);			
								vez.addEdgeterm(ii,jj,kk,ii,jj+1,kk,tempd*-2);
							}
//printf("reach\n");
							if(kk<voxelsize-1&&k+kk<_layers-1)
							{
								double tempd = cost(im(i+ii,j+jj,k+kk),im(i+ii,j+jj,k+kk+1));

								vex.addDataterm(ii,jj,kk,tempd);
								vex.addDataterm(ii,jj,kk+1,tempd);			
								vex.addEdgeterm(ii,jj,kk,ii,jj,kk+1,tempd*-2);

								vey.addDataterm(ii,jj,kk,tempd);
								vey.addDataterm(ii,jj,kk+1,tempd);			
								vey.addEdgeterm(ii,jj,kk,ii,jj,kk+1,tempd*-2);

								vez.addDataterm(ii,jj,kk,tempd);
								vez.addDataterm(ii,jj,kk+1,tempd);			
								vez.addEdgeterm(ii,jj,kk,ii,jj,kk+1,tempd*-2);
							}

							if(ii==0&&i+ii>0)
							{
								vex.addDataterm(ii,jj,kk,xsign*cost(im(i+ii,j+jj,k+kk),im(i+ii-1,j+jj,k+kk)));
							}
							else if(ii==voxelsize-1&&i+ii<_width-1)
							{
								vex.addDataterm(ii,jj,kk,-xsign*cost(im(i+ii,j+jj,k+kk),im(i+ii+1,j+jj,k+kk)));
							}

							if(jj==0&&j+jj>0)
							{
								vey.addDataterm(ii,jj,kk,ysign*cost(im(i+ii,j+jj,k+kk),im(i+ii,j+jj-1,k+kk)));
							}
							else if(jj==voxelsize-1&&j+jj<_height-1)
							{
								vey.addDataterm(ii,jj,kk,-ysign*cost(im(i+ii,j+jj,k+kk),im(i+ii,j+jj+1,k+kk)));
							}

							if(kk==0&&k+kk>0)
							{
								vez.addDataterm(ii,jj,kk,zsign*cost(im(i+ii,j+jj,k+kk),im(i+ii,j+jj,k+kk-1)));
							}
							else if(kk==voxelsize-1&&k+kk<_layers-1)
							{
								vez.addDataterm(ii,jj,kk,-zsign*cost(im(i+ii,j+jj,k+kk),im(i+ii,j+jj,k+kk+1)));
							}

						}
					}
				}

				//getchar();
				//printf("done0\n");
				vex.minimize();
				vey.minimize();
				vez.minimize();
				//getchar();
				//printf("done0\n");
				for(int ii=0;ii<voxelsize&&i+ii<_width;ii++)
				{
//printf("in\n");
					for(int jj=0;jj<voxelsize&&j+jj<_height;jj++)
					{
						for(int kk=0;kk<voxelsize&&k+kk<_layers;kk++)
						{
//printf("iijjkk=%d %d %d %d\n",ii,jj,kk,vex.getLabel(ii,jj,kk));
							xlabel[ii][jj][kk] = (2*vex.getLabel(ii,jj,kk)-1)*(i/voxelsize*2+1)+xsign;
							ylabel[ii][jj][kk] = (2*vey.getLabel(ii,jj,kk)-1)*(j/voxelsize*2+1)+ysign;
							zlabel[ii][jj][kk] = (2*vez.getLabel(ii,jj,kk)-1)*(k/voxelsize*2+1)+zsign;
							_label[i+ii][j+jj][k+kk] = xlabel[ii][jj][kk]*8+ylabel[ii][jj][kk]*4;//+zlabel[i+ii][j+jj][k+kk];
							//_label[i+ii][j+jj][k+kk] = zlabel[i+ii][j+jj][k+kk];

						}
					}
				}
				
				vex.clear();
				vey.clear();
				vez.clear();

				//printf("done2\n");//exit(0);

			

			}
		}
//printf("%d %d %d\n",i,xsign,i+xsign);getchar();

		//if(i>200)
		//{
		//	getchar();

		//	for(int i=0;i<_width;i++)
		//	{
		//		for(int j=0;j<_height;j++)
		//		{
		//			_label[i][j].clear();
		//		}
		//		_label[i].clear();
		//	}
		//	_label.clear();

		//	getchar();
		//	exit(0);
		//}
	}


	return;
}

template<class T>
int gridSV<T>::getLabel(int i, int j, int k)
{
	return _label[i][j][k];
}

#endif