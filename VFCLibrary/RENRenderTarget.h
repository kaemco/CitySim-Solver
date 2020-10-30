#ifndef _INC_RENRENDERTARGET_INCLUDED
#define _INC_RENRENDERTARGET_INCLUDED

#include <memory>
#include <iostream>
#include <cmath>

//#define round(x) (x<0?std::ceil((x)-0.5):std::floor((x)+0.5))

struct SViewPort
{
	SViewPort() : width(0), height(0) {}
    unsigned int width;
    unsigned int height;

	bool IsInside(const int queryX, const int queryY) const {
		return ( (queryX>=0) && (queryX<=static_cast<int>(width)-1) &&
				(queryY>=0) && (queryY<=static_cast<int>(height)-1) );
	}
};

/*
 * This class implements a render target with a colour and depth buffer.  Note that this class
 * does not actually do anything with the depth buffer - it provides it so that
 * some external rasteriser can make use of it if desired.
 */
class RENRenderTarget
{
public:
	RENRenderTarget(void);
	virtual ~RENRenderTarget(void);

	/// Return pointer to the colour buffer
	unsigned int* GetColourBuffer();

    void SetPixels(float *pixels, unsigned int maxSurfaceId) {
        float minPixel = 1e9;
        float maxPixel = 0;
        for (unsigned int i=0; i<GetWidth()*GetHeight(); i++) {
            m_colourBuffer[i]=static_cast<unsigned int>(round(pixels[i]*static_cast<float>(maxSurfaceId)));
            pixels[i]=round(pixels[i]*static_cast<float>(maxSurfaceId));
            if (pixels[i] < minPixel) minPixel = pixels[i];
            if (pixels[i] > maxPixel) maxPixel = pixels[i];
        }
        std::cerr << "Min colour: " << minPixel << "\tMax colour: " << maxPixel << "\tmaxSurfaceColour: " << maxSurfaceId << std::endl;
    }

	void SetPixel(unsigned int x, unsigned int y, unsigned int val) {
		m_colourBuffer[XYtoIndex(x,y)]=val;
	}
	void SetPixel(unsigned int x, unsigned int y, unsigned int val, float depth) {
		m_colourBuffer[XYtoIndex(x,y)]=val;
		m_depthBuffer[XYtoIndex(x,y)]=depth;
	}
	float GetDepth(unsigned int x, unsigned int y) {
		return m_depthBuffer[XYtoIndex(x,y)];
	}
	unsigned int GetPixel(unsigned int x, unsigned int y) {
		return m_colourBuffer[XYtoIndex(x,y)];
	}

	/// Return pointer to the depth buffer
	float* GetDepthBuffer();

	/// Return the 'stride' i.e. number of pixels between adjacent rows
	int GetStride() const;

	/// The width of the render area
	unsigned int GetWidth() const;

	/// The height of the render area
	unsigned int GetHeight() const;

	/// Clear the colour buffer to the specified colour
	void ClearColourBuffer(unsigned int colour);
	void ClearColourBuffer();

	/// Clear the depth buffer to the specified depth
	void ClearDepthBuffer(float depth = 1);

	/// Set the size of the render target
	void SetSize(unsigned int width, unsigned int height);

	const SViewPort& GetViewPort() const { return m_viewport; }

private:
	unsigned int XYtoIndex(unsigned int x, unsigned int y) {
		return x+y*GetStride();
	}

	/// \todo - can we make an equivalent auto_ptr for these?
	unsigned int* m_colourBuffer;
	float* m_depthBuffer;

	SViewPort m_viewport;
};

#endif //_INC_RENRENDERTARGET_INCLUDED
