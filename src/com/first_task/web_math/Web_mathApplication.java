package com.first_task.web_math;

import com.gui.GraphicsDraw;
import com.utils.math.Complex;
import com.vaadin.Application;
import com.vaadin.data.validator.RegexpValidator;
import com.vaadin.ui.Button;
import com.vaadin.ui.Button.ClickEvent;
import com.vaadin.ui.Button.ClickListener;
import com.vaadin.ui.Form;
import com.vaadin.ui.Label;
import com.vaadin.ui.MenuBar;
import com.vaadin.ui.TextField;
import com.vaadin.ui.VerticalLayout;
import com.vaadin.ui.Window;

public class Web_mathApplication extends Application 
{
	
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
		menu.addItem("Main",null);
		menu.addItem("Math",null);
		menu.setSizeFull();
		fld.setValue("1");
		line.setImmediate(true);
		btn.addListener(new ClickListener() 
		{
			
			@Override
			public void buttonClick(ClickEvent event) 
			{
				line.setValue(Double.valueOf(fld.getValue().toString()));
				line.drawSin();
				line.requestRepaint();
				testForm.requestRepaint();
			}
		});
		
		
		
		
		fld.addValidator(new RegexpValidator(Complex.getRegexp(),true, "Eneter complex value"));
		fld.setImmediate(true);
		
		Window mainWindow = new Window("Web_math Application");
		mainWindow.setCaption("Web Math Application");
		mainWindow.addComponent(menu);
		
		
		mainWindow.addComponent(testForm);
		VerticalLayout vl = new VerticalLayout();
		vl.setCaption("vl");
		testForm.setLayout(vl);
		vl.addComponent(label);
		testForm.addField("InputComplex", fld);
		vl.addComponent(btn);
		vl.addComponent(line);


		setMainWindow(mainWindow);
	}

}
