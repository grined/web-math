package com.gui;

import java.util.ArrayList;
import java.util.List;

import com.bibounde.vprotovis.LineChartComponent;
import com.bibounde.vprotovis.common.Point;

public class GraphicsDraw extends  LineChartComponent
{
	private double k=1;
	public GraphicsDraw()
	{
		super();
		this.addSerie("",new Point[] {new Point(-1d,-1d),new Point(2d,2d)});
		this.setChartWidth(600d);
		this.setChartHeight(400d); 

		this.setMarginLeft(20d);

		this.setMarginTop(20d);

		this.setMarginBottom(20d);

        this.setXAxisVisible(true);
        this.setXAxisLabelVisible(true);
        this.setXAxisLabelStep(0.5d);
        this.setXAxisGridVisible(true);
       
        
        
        this.setYAxisVisible(true);
        this.setYAxisLabelVisible(true);
        this.setYAxisLabelStep(0.5d);
        this.setYAxisGridVisible(true);

        this.setTooltipEnabled(false);
	}
	
	private List<Point> calculateSin() 
	{
		ArrayList<Point> points = new ArrayList<>();
		double x = -3.5d;
		while (x < 3.5d)
		{
			points.add(new Point(x,k*Math.sin(x)));
			x+=0.01;
		}
		return points;
	}
	
	public void drawSin() 
	{
		this.clearSeries();
		List points = calculateSin();
		this.addSerie("sin(x)",(Point[]) points.toArray(new Point[points.size()]));
	}
	
	public void setValue(double k)
	{
		this.k = k;
	}
	
}
