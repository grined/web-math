package com.first_task.web_math;

import com.gui.GraphicsDraw;
import com.vaadin.Application;
import com.vaadin.ui.AbsoluteLayout;
import com.vaadin.ui.Button;
import com.vaadin.ui.Form;
import com.vaadin.ui.Label;
import com.vaadin.ui.MenuBar;
import com.vaadin.ui.Panel;
import com.vaadin.ui.TextField;
import com.vaadin.ui.Window;

public class Web_mathApplication extends Application 
{
	private static String ERR_COMPLEX = "Enter complex value";
	Panel mainPanel = new Panel();
	
	Button naviMain = 		new Button("Main");
	Button naviExample = 	new Button("Example");
	Button naviContacts = 	new Button("Contacts");
	
	TextField tfE1 = new TextField("E1");
	TextField tfE2 = new TextField("E2");
	TextField tfE3 = new TextField("E3");

	TextField tfnu1 = new TextField("nu1");
	TextField tfnu2 = new TextField("nu2");
	TextField tfnu3 = new TextField("nu3");

	TextField tfa1 = new TextField("a1");
	TextField tfa2 = new TextField("a2");
	TextField tfa3 = new TextField("a3");
	
	TextField tfK1 = new TextField("K1");
	TextField tfK2 = new TextField("K2");
	TextField tfK3 = new TextField("K3");
	
	TextField tfk= new TextField("k");
	TextField tfhp= new TextField("hp");
	TextField tfn= new TextField("n");
	TextField tfT01= new TextField("T01");
	TextField tfT02= new TextField("T02");
	
	TextField tfH= new TextField("H (use ';' )");
	TextField tfh= new TextField("h (use ';' )");
	
	
	
	
	Button btn = new Button("Press me");
	Label label = new Label("Hello Vaadin user");
	MenuBar menu = new MenuBar();
	GraphicsDraw line = new GraphicsDraw();
	TextField fld = new TextField("Enter value");
	Form testForm = new Form();
	
	@SuppressWarnings("serial")
	@Override
	public void init() 
	{
		setTheme("math_cover_theme");
        
 
        
//		fld.setValue("1");
//		line.setImmediate(true);
//		btn.addListener(new ClickListener() 
//		{
//			
//			@Override
//			public void buttonClick(ClickEvent event) 
//			{
//				line.setValue(Double.valueOf(fld.getValue().toString()));
//				line.drawSin();
//				line.requestRepaint();
//				testForm.requestRepaint();
//			}
//		});

		
//TODO: for feautyre
//		fld.addValidator(new RegexpValidator(Complex.getRegexp(),true, "Eneter complex value"));
//		fld.setImmediate(true);
		
		Window mainWindow = new Window("Web_math Application");
		mainWindow.setCaption("Web Math Application");
		setMainWindow(mainWindow);
		
		mainWindow.addComponent(mainPanel);
//		mainWindow.setSizeUndefined();
		mainPanel.setScrollable(true);
//		mainPanel.setSizeFull();
		
		AbsoluteLayout mainLayout = new AbsoluteLayout();
		mainLayout.setSizeFull();
		mainLayout.setWidth("100%");
		mainLayout.setHeight("100%");
				
		mainPanel.addComponent(mainLayout);
		// Navi
		AbsoluteLayout layNavi = new AbsoluteLayout();
		mainLayout.addComponent(layNavi, xyCSS(200,100));

		int x = 230,y = 20;
		naviMain.setWidth("80px");
		naviExample.setWidth("80px");
		naviContacts.setWidth("80px");
		
		layNavi.addComponent(naviMain, xyCSS(x,y));
		layNavi.addComponent(naviExample, xyCSS(x+150,y));
		layNavi.addComponent(naviContacts, xyCSS(x+300,y));
		
		// MathLayout
		AbsoluteLayout layMath = new AbsoluteLayout();
		mainLayout.addComponent(layMath, xyCSS(200,200));//TODO:SWITCH

		int xm = 20, ym=20;
		layMath.addComponent(new Label("Enter values:"),xyCSS(xm,ym));
		
		layMath.addComponent(tfE1,xyCSS(xm,ym+30));
		layMath.addComponent(tfE2,xyCSS(xm+100,ym+30));
		layMath.addComponent(tfE3,xyCSS(xm+200,ym+30));
		
//		for ( it: layMath.getComponentIterator())
//		{
//			
//		}
		tfE1.setWidth("80px");
		tfE2.setWidth("80px");
		tfE3.setWidth("80px");
		
//		vl.addComponent(label);
//		vl.addComponent(fld);
//		//testForm.addField("InputComplex", fld);
//		vl.addComponent(btn);
//		vl.addComponent(line);
//		
//		vl.setComponentAlignment(label, 100, 50);
//		vl.setComponentAlignment(btn, 100, 50);
//		vl.setComponentAlignment(fld, 100, 50);
//		vl.setComponentAlignment(line, 100, 50);

		setMainWindow(mainWindow);
	}
	
	public String xyCSS(int x, int y)
	{
		String xs = String.valueOf(x);
		String ys = String.valueOf(y);
		
		return "top:"+ys+".0px;left:"+xs+".0px;";
	}

}
