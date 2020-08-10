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
package quantumVITAS;

import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Order;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.TestMethodOrder;
import org.junit.jupiter.api.MethodOrderer;
import org.testfx.api.FxRobotException;
import org.testfx.service.query.EmptyNodeQueryException;
import org.testfx.util.WaitForAsyncUtils;
import javafx.collections.ObservableList;
import javafx.scene.Node;
import javafx.scene.Parent;
import javafx.scene.control.ComboBox;
import javafx.scene.control.RadioButton;
import javafx.scene.control.Tab;
import javafx.scene.control.TabPane;
import javafx.scene.control.TextField;
import javafx.scene.control.TreeItem;
import javafx.scene.control.TreeTableView;
import javafx.scene.input.KeyCode;
import project.ProjectCalcLog;

@TestMethodOrder(MethodOrderer.OrderAnnotation.class)
public class FirstMainWindowTest extends MainWindowTest{
	
	@Test
	@Order(1)
	public void testButtonException() {
		Assertions.assertThrows(FxRobotException.class, () -> {
		    clickOn("#anyButtonNonExisting");
		  });
		Assertions.assertThrows(EmptyNodeQueryException.class, () -> {
			lookup("#anyLabelNonExisting").queryLabeled();
		  });
	}
	
	@Test
	@Order(2)
	public void testSetupWorkSpace() throws TimeoutException {
		clickOn("#buttonOpenWorkSpace");
		TextField textWorkSpace = (TextField) lookup("#textWorkSpace").query();
		Assertions.assertTrue(textWorkSpace.getText().contains("testfx"),"textWorkSpace string should contain 'testfx' in the folder name");
	}
	
	@Test
	@Order(3)
	public void testCreateProject() throws TimeoutException {
		createProject();
	}
	
	
	@Test
	@Order(4)
	public void testAddCalculation() {
		TreeTableView<ProjectCalcLog> projectTree = lookup("#projectTree").query();
		TabPane tabPaneRight = (TabPane) lookup("#idtabPaneRight").queryParent();
		ObservableList<Tab> tabList;
		
		clickOn("#calcMain");
		clickOn("#calcScf");
		
		Assertions.assertTrue(containsRow(projectTree,"testProject0","SCF_1"));
		tabList = tabPaneRight.getTabs();
		Assertions.assertTrue(tabList.size()==2);//geo and scf tab
		//Assertions.assertTrue(comboCalculation.getSelectionModel().getSelectedItem().toLowerCase().contains("scf"));
		
		clickOn("#calcMain");
		clickOn("#calcOpt");
		Assertions.assertTrue(containsRow(projectTree,"testProject0","OPT_1"));
		tabList = tabPaneRight.getTabs();
		Assertions.assertTrue(tabList.size()==3);//geo, scf and opt

		clickOn("#calcMain");
		clickOn("#calcMd");
		Assertions.assertTrue(containsRow(projectTree,"testProject0","BOMD_1"));
		tabList = tabPaneRight.getTabs();
		Assertions.assertTrue(tabList.size()==3);//geo, scf and md
	}
	
	@Test
	@Order(5)
	public void testToggleToGeo() {
		clickOn("Geometry"); 
		
		TabPane tabPaneRight = (TabPane) lookup("#idtabPaneRight").queryParent();
		ObservableList<Tab> tabList = tabPaneRight.getTabs();
		Assertions.assertTrue(tabList.size()==1);//only geo tab
		Assertions.assertTrue(tabList.get(0).getText().toLowerCase().contains("geo"));
	}
	
	@Test
	@Order(6)
	public void testGeoBravais() {
		Node nd = lookup("#titlePaneCell").query();
		Node ndNew = from(nd).lookup(".title").query();
		clickOn(ndNew);
		
		//sleep(1000);
		ComboBox<String> ibravCombo = lookup("#ibravCombo").queryComboBox();
		selectComboBox(ibravCombo, 2);
		
		TextField aField = (TextField) lookup("#aField").queryTextInputControl();
		clickOn(aField).type(KeyCode.DIGIT5,KeyCode.DECIMAL,KeyCode.DIGIT4);
		//sleep(1000);
	}
	
	@Test
	@Order(7)
	public void testGeoAtoms() {
		Node nd = lookup("#titlePaneCell").query();
		Node ndNew = from(nd).lookup(".title").query();
		clickOn(ndNew);
		
		nd = lookup("#titlePaneAtoms").query();
		ndNew = from(nd).lookup(".title").query();
		clickOn(ndNew);
		
		sleep(1000);
		TextField textElem = (TextField) from(nd).lookup("#textElem").queryTextInputControl();
		clickOn(textElem).type(KeyCode.CAPS,KeyCode.S,KeyCode.CAPS,KeyCode.I);//Si
		TextField textX = (TextField) from(nd).lookup("#textX").queryTextInputControl();
		clickOn(textX).type(KeyCode.DIGIT0);//0
		TextField textY = (TextField) from(nd).lookup("#textY").queryTextInputControl();
		clickOn(textY).type(KeyCode.DIGIT0);//0
		TextField textZ = (TextField) from(nd).lookup("#textZ").queryTextInputControl();
		clickOn(textZ).type(KeyCode.DIGIT0);//0.0

		Node addButton = from(nd).lookup("#buttonAddEnd").query();
		Node clearButton = from(nd).lookup("#clearInput").query();
		
		clickOn(addButton);
		clickOn(clearButton);
		
		clickOn(textElem).type(KeyCode.CAPS,KeyCode.S,KeyCode.CAPS,KeyCode.I);//Si
		clickOn(textX).type(KeyCode.DIGIT0,KeyCode.DECIMAL,KeyCode.DIGIT2,KeyCode.DIGIT5);//0.25
		clickOn(textY).type(KeyCode.DIGIT0,KeyCode.DECIMAL,KeyCode.DIGIT2,KeyCode.DIGIT5);//0.25
		clickOn(textZ).type(KeyCode.DIGIT0,KeyCode.DECIMAL,KeyCode.DIGIT2,KeyCode.DIGIT5);//0.25
		
		clickOn(addButton);
		
		Node root3DVbox = lookup("#root3DVbox").query();
		RadioButton radio3Cart = (RadioButton) from(root3DVbox).lookup("#radio3Cart").query();
		clickRadioButton(radio3Cart);
		
		TextField tfx = (TextField) from(root3DVbox).lookup("#tfx").queryTextInputControl();
		TextField tfy = (TextField) from(root3DVbox).lookup("#tfy").queryTextInputControl();
		TextField tfz = (TextField) from(root3DVbox).lookup("#tfz").queryTextInputControl();
		clickOn(tfx).type(KeyCode.DIGIT1);//1
		clickOn(tfy).type(KeyCode.DIGIT1);//1
		clickOn(tfz).type(KeyCode.DIGIT1);//1
		
		Node btUpd = from(root3DVbox).lookup("#btUpd").query();//update button
		clickOn(btUpd);
		

		RadioButton radio2Cryst = (RadioButton) from(root3DVbox).lookup("#radio2Cryst").query();
		clickRadioButton(radio2Cryst);
		clickOn(tfx).type(KeyCode.BACK_SPACE,KeyCode.DIGIT2);//2
		clickOn(tfy).type(KeyCode.BACK_SPACE,KeyCode.DIGIT2);//2
		clickOn(tfz).type(KeyCode.BACK_SPACE,KeyCode.DIGIT2);//2
		clickOn(btUpd);
		
	}
	
	@Test
	@Order(8)
	public void testGoToCalculation() {
		Node treeNode = lookup("#projectTree").query();
		TabPane tabPaneRight = (TabPane) lookup("#idtabPaneRight").queryParent();
		
		
		Node nodeScf1 = from(treeNode).lookup("SCF_1").query();
		//sleep(1000);
		clickOn(nodeScf1);
		checkRightPaneTabs("scf", 0, tabPaneRight);
		

		nodeScf1 = from(treeNode).lookup("OPT_1").query();
		clickOn(nodeScf1);
		checkRightPaneTabs("opt", 1, tabPaneRight);
		

		nodeScf1 = from(treeNode).lookup("BOMD_1").query();
		clickOn(nodeScf1);
		checkRightPaneTabs("md", 2, tabPaneRight);

	}
	
//	@Test
//	@Order(9)
//	public void testCheckBoxResetAll() {
//		//VBox vboxScfStandard = (VBox)lookup("#vboxStandard").query();
//		//printQueryAll("#checkResetAll");
//		
//		
//		ComboBox<String> comboCalculation = lookup("#comboCalculation").queryComboBox();
//		TabPane tabPaneRight = (TabPane) lookup("#idtabPaneRight").queryParent();
//		
//		//go to scf
//		selectComboBox(comboCalculation, 0);
//		
//		Tab tabSelected = tabPaneRight.getSelectionModel().getSelectedItem();
//		
//		Node nodeCheckResetAll = from(tabSelected.getContent()).lookup("#checkResetAll").query();
//		//printQueryAll(tabSelected.getContent(),"#checkResetAll");
//		
//		CheckBox cb = (CheckBox) nodeCheckResetAll;
//		Node nodeScfStandard = from(tabSelected.getContent()).lookup("#standardPane").query();
//		//sleep(1000);
//		clickOn(nodeScfStandard);
//		//sleep(1000);
//		//clickOn(nodeCheckResetAll);
//		clickCheckBox(cb);
//		//sleep(1000);
//		Assertions.assertTrue(cb.isSelected(),"Reset all checkbox should be selected!");
//	}
//	
//	
//	
//	@Test
//	@Order(10)
//	public void testSaveProjectButton() {
//		ComboBox<String> comboProject = lookup("#comboProject").queryComboBox();
//		String projectFromCombo = comboProject.getSelectionModel().getSelectedItem();
//		
//		TabPane tabPaneRight = (TabPane) lookup("#workSpaceTabPane").queryParent();
//		String projectFromTab = tabPaneRight.getSelectionModel().getSelectedItem().getText();
//		
//		Assertions.assertEquals(projectFromCombo, projectFromTab,"|"+projectFromCombo+"|"+projectFromTab+"|");
//
//		clickOn("#saveProjectButton");
//	}
//	
//	@Test
//	@Order(11)
//	public void testOpenCloseProject() {
//		ComboBox<String> comboCalculation = lookup("#comboCalculation").queryComboBox();
//		ComboBox<String> comboProject = lookup("#comboProject").queryComboBox();
//		String projectFromCombo = comboProject.getSelectionModel().getSelectedItem();
//		
//		//printQueryAll("#workSpaceTabPane");
//		
//		TabPane workSpaceTabPane = (TabPane) lookup("#workSpaceTabPane").query();
//		
//		Node mainScrollLeft = lookup("#mainScrollLeft").query();
//		Node projectTreeItem = from(mainScrollLeft).lookup(projectFromCombo).query();
//		Node buttonCloseSelected = from(mainScrollLeft).lookup("#buttonCloseSelected").query();
//		Node buttonOpenSelected = from(mainScrollLeft).lookup("#buttonOpenSelected").query();
//		//printQueryAll(mainScrollLeft,projectFromCombo);
//		
//		//close the project
//		clickOn(projectTreeItem);
//		clickOn(buttonCloseSelected);
//		
//		//sleep(1000);
//		//now there should still be the treeitem, but there should not be any tab
//		projectTreeItem = from(mainScrollLeft).lookup(projectFromCombo).query();
//		int numTabs = workSpaceTabPane.getTabs().size();
//		Assertions.assertTrue(numTabs==0,projectFromCombo+Integer.toString(numTabs));
//		Assertions.assertTrue(comboCalculation.getItems().isEmpty());
//		
//		//now open the project again
//		clickOn(projectTreeItem);
//		clickOn(buttonOpenSelected);
//		
//		Tab selectedTab = workSpaceTabPane.getSelectionModel().getSelectedItem();
//		boolean fl = (selectedTab!=null) && (projectFromCombo.equals(selectedTab.getText()));
//		Assertions.assertTrue(fl);
//		Assertions.assertTrue(!comboCalculation.getItems().isEmpty());
//	}
//	
//	@Test
//	@Order(12)
//	public void testCreateProjectMore() throws TimeoutException {
//		createProject();
//		createProject();
//	}
	
	private void printQueryAll(String queryStr) {
		Set<Node> allNodes = lookup(queryStr).queryAll();
		String msg="\nNodes\n";
		for(Node nd : allNodes) {
			msg+=nd.toString()+nd.isVisible()+",\n";
			for(Node nd_child : ((Parent)nd).getChildrenUnmodifiable()) {
				msg+="---Children:"+nd_child.toString()+","+nd_child.getStyleClass()+",\n";
			}
		}
		msg+="Parents:\n";
		for(Node nd : allNodes) {
			msg+=nd.getParent().toString()+",\n";
		}
		Assertions.assertTrue(false,msg);
	}
	private void printQueryAll(Node rootNode, String queryStr) {
		Set<Node> allNodes = from(rootNode).lookup(queryStr).queryAll();
		String msg="\nNodes\n";
		for(Node nd : allNodes) {
			msg+=nd.toString()+nd.isVisible()+",\n";
			for(Node nd_child : ((Parent)nd).getChildrenUnmodifiable()) {
				msg+="---Children:"+nd_child.toString()+","+nd_child.getStyleClass()+",\n";
			}
		}
		msg+="Parents:\n";
		for(Node nd : allNodes) {
			msg+=nd.getParent().toString()+",\n";
		}
		Assertions.assertTrue(false,msg);
	}
	private void checkRightPaneTabs(String selectItem, int iExpected, TabPane tabPaneRight) {

		ObservableList<Tab> tabList = tabPaneRight.getTabs();
		
		if(selectItem.contains("scf")&&!selectItem.contains("nscf")) {Assertions.assertTrue(tabList.size()==2,"scf"+Integer.toString(tabList.size()));}
		if(selectItem.contains("opt")) {Assertions.assertTrue(tabList.size()==3,"opt"+Integer.toString(tabList.size()));}
		if(selectItem.contains("md")) {Assertions.assertTrue(tabList.size()==3,"md"+Integer.toString(tabList.size()));}
	}
	private <T> void selectComboBox(ComboBox<T> cb, int i) {
		if(i<0 || i>cb.getItems().size()) {Assertions.assertTrue(false,"Index out of bound");return;}
		
		Node nd_tmp=null;
		String msg="";
		int maxTry = 5;
		int tryCount=0;
		do {
			tryCount++;
			clickOn(cb);
			
			Set<Node> allListCells = lookup(cb.getItems().get(i).toString()).queryAll();
					
			msg="";
			
			for(Node nd : allListCells) {
				//if(nd.getStyleClass().contains("list-cell")) {
				if(nd.getStyleClass().contains("list-cell")) {
					nd_tmp = nd;
					msg+=nd.toString()+",yes,\n";
					
				}
				else {
					msg+=nd.toString()+","+nd.getParent().toString()+","+nd.getParent().getStyleClass()+",no,\n";
				}
			}
		}while(nd_tmp==null && tryCount<maxTry);
		
		Assertions.assertTrue(nd_tmp!=null,cb.getItems().get(i).toString()+"\n"+msg);
		clickOn(nd_tmp);
	}
	private Node lookupById(final String controlId) throws Exception {
		WaitForAsyncUtils.waitFor(2, TimeUnit.SECONDS, new Callable<Boolean>() {
		    @Override
		    public Boolean call() throws Exception {
		        return lookup(controlId).query() != null;
		    }
		});
		return lookup(controlId).query();
	}
	private void createProject() throws TimeoutException {
		clickOn("#createProject");
		//Button createProject = lookup("#createProject").queryButton();
		TreeTableView<ProjectCalcLog> projectTree = lookup("#projectTree").queryAs(TreeTableView.class);
		ObservableList<TreeItem<ProjectCalcLog>> obs = projectTree.getRoot().getChildren();
		for(int i=0;i<obs.size();i++) {
			Assertions.assertTrue(obs.get(i).getValue().getProject().contains("testProject"),
					"projectTree selected item should contain 'testProject'");
		}
	}
//	private void clickCheckBox(CheckBox cb) {
//		for (Node child : cb.getChildrenUnmodifiable()) {
//	        if (child.getStyleClass().contains("box")) {
//	            clickOn(child);
//	        }
//	    }
//	}
	private void clickRadioButton(RadioButton rb) {
		for (Node child : rb.getChildrenUnmodifiable()) {
	        if (child.getStyleClass().contains("radio")) {
	            clickOn(child);
	        }
	    }
	}
	private boolean containsRow(TreeTableView<ProjectCalcLog> projectTree,String strProj, String strCalc) {
		TreeItem<ProjectCalcLog> rootTree = projectTree.getRoot();
		for(TreeItem<ProjectCalcLog> tiProj: rootTree.getChildren()) {
			String strTmpProj = tiProj.getValue().getProject();
			for(TreeItem<ProjectCalcLog> tiCalc: tiProj.getChildren()) {
				if(strCalc.equals(tiCalc.getValue().getCalculation()) && strProj.equals(strTmpProj)) {return true;}
			}
		}
		return false;
	}
	
}
