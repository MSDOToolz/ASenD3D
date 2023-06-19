#ifndef FACE
#define FACE

class Face {
	private:
	    int numNds;
	    int *locNodes;
		int *globNodes;
		bool onSurf;
		Face *next;
		
	public:
	    Face(int newNumNds);
		
		void setNode(int place, int locNd, int globNd);
		
		void setOnSurface(bool newOnSurf);
		
		void setNext(Face *newNext);
		
		void sortedNodes(int srtNds[]);
		
		int getNumNds();
		
		bool onSurface();
		
		int getLowNd();
		
		Face* getNext();
		
		void destroy();
};

class FaceList {
	private:
	    Face *first;
		Face *last;
		int len;
		
	public:
	    FaceList();
		
		void addFace(Face *newFc);
		
		void addIfAbsent(Face *newFc);
		
		Face* getFirst();
		
		void destroy();
};

#endif