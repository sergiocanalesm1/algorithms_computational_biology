package org.myorg.mycode;

public class Atom {
	private String atomName;
	private String atomType;
	private String resName;
	private int resNum;
	private double x_cord;
	private double y_cord;
	private double z_cord;
	private int  atomId;
	
	public Atom(String atomName, String atomType,String resName,	 int resNum,	 double x_cord,	 double y_cord,	 double z_cord,	 int  atomId) {
		this.setAtomName(atomName);
		this.setAtomType(atomType);
		this.setResName(resName);
		this.setResNum(resNum);
		this.setX_cord(x_cord);
		this.setY_cord(y_cord);
		this.setZ_cord(z_cord);
		this.setAtomId(atomId);	
	}



	public String getAtomType() {
		return atomType;
	}

	public void setAtomType(String atomType) {
		this.atomType = atomType;
	}

	/**
	 * @return the resName
	 */
	public String getResName() {
		return resName;
	}

	/**
	 * @param resName the resName to set
	 */
	public void setResName(String resName) {
		this.resName = resName;
	}

	/**
	 * @return the resNum
	 */
	public int getResNum() {
		return resNum;
	}

	/**
	 * @param resNum the resNum to set
	 */
	public void setResNum(int resNum) {
		this.resNum = resNum;
	}

	/**
	 * @return the x_cord
	 */
	public double getX_cord() {
		return x_cord;
	}

	/**
	 * @param x_cord the x_cord to set
	 */
	public void setX_cord(double x_cord) {
		this.x_cord = x_cord;
	}

	/**
	 * @return the y_cord
	 */
	public double getY_cord() {
		return y_cord;
	}

	/**
	 * @param y_cord the y_cord to set
	 */
	public void setY_cord(double y_cord) {
		this.y_cord = y_cord;
	}

	/**
	 * @return the z_cord
	 */
	public double getZ_cord() {
		return z_cord;
	}

	/**
	 * @param z_cord the z_cord to set
	 */
	public void setZ_cord(double z_cord) {
		this.z_cord = z_cord;
	}

	/**
	 * @return the atomId
	 */
	public int getAtomId() {
		return atomId;
	}

	/**
	 * @param atomId the atomId to set
	 */
	public void setAtomId(int atomId) {
		this.atomId = atomId;
	}



	/**
	 * @return the atomName
	 */
	public String getAtomName() {
		return atomName;
	}



	/**
	 * @param atomName the atomName to set
	 */
	public void setAtomName(String atomName) {
		this.atomName = atomName;
	}
	
	
	
	

}
