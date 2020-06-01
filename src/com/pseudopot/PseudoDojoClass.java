/*******************************************************************************
 * Copyright (c) 2020 Haonan Huang.
 *
 *     This file is part of QuantumVITAS (Quantum Visualization Interactive 
 *     Toolkit for Ab-initio Simulations).
 *
 *     QuantumVITAS is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or any 
 *     later version.
 *
 *     QuantumVITAS is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with QuantumVITAS.  If not, see <https://www.gnu.org/licenses/gpl-3.0.txt>.
 *******************************************************************************/

package com.pseudopot;

import java.io.File;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;

import com.consts.Constants.EnumFunctional;
import com.pseudopot.PseudoDojoEnum.*;

import javafx.scene.control.Alert;
import javafx.scene.control.Alert.AlertType;

public class PseudoDojoClass extends PseudoPotential{

	private EnumFunctional typeFunctional;
	private String precString;
	private boolean isRelativ;//whether or not fully relativistic
	
	public PseudoDojoClass() {
		super(EnumPseudoPotLib.PSEUDODOJO, true);
		precisionList.add("Standard");precisionList.add("Stringent");
		functionalList.add(EnumFunctional.PBE);
		functionalList.add(EnumFunctional.PBESOL);
		functionalList.add(EnumFunctional.LDA);
		
	}
	public String getFile(String element) { //element name in the format of e.g. "He" 
		if (typeFunctional==null || !functionalList.contains(typeFunctional)) return null;
		if (precString==null || !precisionList.contains(precString)) return null;
		
		String st1 =  getString(element, "getFolderName");
		String st2 =  getString(element, "getFileName");
		
		if(st1!=null && st2!=null) {
			return st1+File.separator+st2;
		}
		else {return null;}
	}
	public String getString(String element, String methodName) { //element name in the format of e.g. "He" 
		//use reflection to simplify code
		
		if (typeFunctional==null || !functionalList.contains(typeFunctional)) return null;
		if (precString==null || !precisionList.contains(precString)) return null;
		
		Method method;
		try {
			switch(typeFunctional) {
				case PBE:
					if (precString.equals("Standard")) {
						if (isRelativ) {
							fr_pbe_standard ef = fr_pbe_standard.valueOf(element);
							method = ef.getClass().getMethod(methodName);
							return (String) method.invoke(ef);}
						else {
							sr_pbe_standard ef = sr_pbe_standard.valueOf(element);
							method = ef.getClass().getMethod(methodName);
							return (String) method.invoke(ef);}
					}
					else if (precString.equals("Stringent")) {
						if (isRelativ) {
							fr_pbe_stringent ef = fr_pbe_stringent.valueOf(element);
							method = ef.getClass().getMethod(methodName);
							return (String) method.invoke(ef);}
						else {
							sr_pbe_stringent ef = sr_pbe_stringent.valueOf(element);
							method = ef.getClass().getMethod(methodName);
							return (String) method.invoke(ef);}
					}
					else return null;
				case PBESOL:
					if (precString.equals("Standard")) {
						if (isRelativ) {
							fr_pbesol_standard ef = fr_pbesol_standard.valueOf(element);
							method = ef.getClass().getMethod(methodName);
							return (String) method.invoke(ef);}
						else {
							sr_pbesol_standard ef = sr_pbesol_standard.valueOf(element);
							method = ef.getClass().getMethod(methodName);
							return (String) method.invoke(ef);}
					}
					else if (precString.equals("Stringent")) {
						if (isRelativ) {
							fr_pbesol_stringent ef = fr_pbesol_stringent.valueOf(element);
							method = ef.getClass().getMethod(methodName);
							return (String) method.invoke(ef);}
						else {
							sr_pbesol_stringent ef = sr_pbesol_stringent.valueOf(element);
							method = ef.getClass().getMethod(methodName);
							return (String) method.invoke(ef);}	
					}
					else return null;
				case LDA:
					if (precString.equals("Standard")) {
						if (isRelativ) {return null;}
						else {
							sr_pw_standard ef = sr_pw_standard.valueOf(element);
							method = ef.getClass().getMethod(methodName);
							return (String) method.invoke(ef);}
					}
					else if (precString.equals("Stringent")) {
						if (isRelativ) {return null;}
						else {
							sr_pw_stringent ef = sr_pw_stringent.valueOf(element);
							method = ef.getClass().getMethod(methodName);
							return (String) method.invoke(ef);}
					}
					else return null;
				default:return null;
			}
		}
		catch (IllegalArgumentException | IllegalAccessException | InvocationTargetException | SecurityException | NoSuchMethodException e) {
			Alert alert = new Alert(AlertType.INFORMATION);
	    	alert.setTitle("Error");
	    	alert.setContentText("Error in PseudoDojoClass.getString()! "+e.getMessage());
	    	alert.showAndWait();
	    	
			return null;
		}
	}
	public EnumFunctional getTypeFunctional() {
		return typeFunctional;
	}

	public void setTypeFunctional(EnumFunctional tf) {
		this.typeFunctional = tf;
	}

	public String getPrecString() {
		return precString;
	}

	public void setPrecString(String tp) {
		//should take null!
		this.precString = tp;
	}

	public boolean isRelativ() {
		return isRelativ;
	}

	public void setRelativ(boolean isRelativ) {
		this.isRelativ = isRelativ;
	}
	
//	public String getFile(String element) { //element name in the format of e.g. "He" 
//	if (typeFunctional==null || !functionalList.contains(typeFunctional)) return null;
//	if (precString==null || !precisionList.contains(precString)) return null;
//	
//	//return getString(element, "getFolderName");
//	
//	try {
//		switch(typeFunctional) {
//			case PBE:
//				if (precString.equals("Standard")) {
//					if (isRelativ) {fr_pbe_standard ef = fr_pbe_standard.valueOf(element);
//					return ef.getFolderName()+File.separator+ef.getFileName();}
//					else {sr_pbe_standard ef = sr_pbe_standard.valueOf(element);
//					return ef.getFolderName()+File.separator+ef.getFileName();}
//				}
//				else if (precString.equals("Stringent")) {
//					if (isRelativ) {fr_pbe_stringent ef = fr_pbe_stringent.valueOf(element);
//					return ef.getFolderName()+File.separator+ef.getFileName();}
//					else {sr_pbe_stringent ef = sr_pbe_stringent.valueOf(element);
//					return ef.getFolderName()+File.separator+ef.getFileName();}
//				}
//				else return null;
//			case PBESOL:
//				if (precString.equals("Standard")) {
//					if (isRelativ) {fr_pbesol_standard ef = fr_pbesol_standard.valueOf(element);
//					return ef.getFolderName()+File.separator+ef.getFileName();}
//					else {sr_pbesol_standard ef = sr_pbesol_standard.valueOf(element);
//					return ef.getFolderName()+File.separator+ef.getFileName();}
//				}
//				else if (precString.equals("Stringent")) {
//					if (isRelativ) {fr_pbesol_stringent ef = fr_pbesol_stringent.valueOf(element);
//					return ef.getFolderName()+File.separator+ef.getFileName();}
//					else {sr_pbesol_stringent ef = sr_pbesol_stringent.valueOf(element);
//					return ef.getFolderName()+File.separator+ef.getFileName();}
//				}
//				else return null;
//			case LDA:
//				if (precString.equals("Standard")) {
//					if (isRelativ) {return null;}
//					else {sr_pw_standard ef = sr_pw_standard.valueOf(element);
//					return ef.getFolderName()+File.separator+ef.getFileName();}
//				}
//				else if (precString.equals("Stringent")) {
//					if (isRelativ) {return null;}
//					else {sr_pw_stringent ef = sr_pw_stringent.valueOf(element);
//					return ef.getFolderName()+File.separator+ef.getFileName();}
//				}
//				else return null;
//			default:return null;
//		}
//	}
//	catch (IllegalArgumentException iae) {
//		return null;
//	}
//
//}

}
