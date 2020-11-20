public class Matrix {
    private static final double EPS = 10e-10;

    private double[][] array;
    private int rows;
    private int columns;

    public Matrix(double[][] array){
        this.array = array;
        this.rows = array.length;
        this.columns = array[0].length;
    }

    //Getters
    public double[][] getArray() {
        return array;
    }
    public double[] getDimension(){
        return new double[]{rows, columns};
    }

    //Calculates algebraic complements
    public Matrix getComplementMatrix() throws IllegalMatrixDimension {
        if(!isSquare()) throw new IllegalMatrixDimension();
        double[][] array = new double[this.rows][this.columns];
        for (int i = 0; i < array.length; i++) {
            for (int j = 0; j < array[0].length; j++) {
                int sign;
                if((i + j) %2 == 0) sign = 1;
                else sign = -1;

                array[i][j] = sign * getMinor(i,j).getDeterminant();
            }
        }
        return new Matrix(array);
    }
    //Calculates determinant of matrix
    public double getDeterminant() throws IllegalMatrixDimension {
        if(!isSquare()) throw new IllegalMatrixDimension();
        if(this.rows < 3) return this.array[0][0] * this.array[1][1] - this.array[0][1] * this.array[1][0];
        double res = 0;
        for (int i = 0; i < this.columns; i++) {
            Matrix minor = this.getMinor(0, i);
            double detMinor = minor.getDeterminant();
            if (i % 2 == 0) res += this.array[0][i] * detMinor;
            else res -= this.array[0][i] * detMinor;
        }
        return res;
    }

    //Actions
    public void multiplyByNumber(double number){
        for (int i = 0; i < this.rows; i++) {
            for (int j = 0; j < this.columns; j++) {
                this.array[i][j] *= number;
            }
        }
        this.complete();
    }
    public void add(Matrix matrix) throws IllegalMatrixDimension {
        if(!canBeAdded(this, matrix)) throw new IllegalMatrixDimension();
        for (int i = 0; i < this.rows; i++) {
            for (int j = 0; j < this.columns; j++) {
                this.array[i][j] += matrix.array[i][j];
            }
        }
        this.complete();
    }
    public void subtract(Matrix matrix) throws IllegalMatrixDimension {
        if(!canBeAdded(this, matrix)) throw new IllegalMatrixDimension();
        for (int i = 0; i < this.rows; i++) {
            for (int j = 0; j < this.columns; j++) {
                this.array[i][j] -= matrix.array[i][j];
            }
        }
        this.complete();
    }
    public void multiplyByMatrix(Matrix matrix) throws IllegalMatrixDimension {
        if(!canBeMultiplied(this, matrix)) throw new IllegalMatrixDimension();

        double[][] array = new double[this.rows][matrix.columns];
        Matrix res = new Matrix(this.rows,matrix.columns);

        for (int i = 0; i < res.rows; i++) {
            for (int j = 0; j < res.columns; j++) {
                for (int k = 0; k < this.columns; k++) {
                    array[i][j] += this.array[i][k] * matrix.array[k][j];
                }
            }
        }

        res.complete(array);
        this.changeTo(res);
    }
    public void divideByMatrix(Matrix matrix) throws IllegalMatrixDimension {
        Matrix invMatrix = matrix.getInverse();
        this.multiplyByMatrix(invMatrix);
    }
    //Group Actions
    public static Matrix summarise(Matrix[] matrices) throws IllegalMatrixDimension {
        if(matrices.length == 0) return null;
        Matrix res = copyMatrix(matrices[0]);
        for (int i = 1; i < matrices.length; i++) {
            res.add(matrices[i]);
        }
        return res;
    }
    public static Matrix multiply(Matrix[] matrices) throws IllegalMatrixDimension {
        if(matrices.length == 0) return null;
        Matrix res = copyMatrix(matrices[0]);
        for (int i = 1; i < matrices.length; i++) {
            res.multiplyByMatrix(matrices[i]);
        }
        return res;
    }
    //Actions for square matrices
    public Matrix getMinor(int row, int column){
        double[][] res = new double[this.rows - 1][this.columns - 1];
        for (int i = 0, k = 0; i < this.rows; i++) {
            if(i == row) continue;
            for (int j = 0, m = 0; j < this.columns; j++) {
                if(j == column) continue;
                res[k][m] = this.array[i][j];
                m++;
            }
            k++;
        }
        return new Matrix(res);
    }
    public Matrix getInverse() throws IllegalMatrixDimension {
        if(!doesInverseExist()) throw new IllegalMatrixDimension();
        double determinant = this.getDeterminant();
        Matrix complementTransposed = this.getComplementMatrix().getTransposed();
        complementTransposed.multiplyByNumber(1/determinant);
        return complementTransposed;
    }
    public Matrix getTransposed(){
        double[][] array = new double[this.columns][this.rows];
        for (int i = 0; i < this.columns; i++) {
            for (int j = 0; j < this.rows; j++) {
                array[i][j] = this.array[j][i];
            }
        }
        return new Matrix(array);
    }
    public Matrix pow(int degree) throws IllegalMatrixDimension {
        if(!isSquare()) throw new IllegalMatrixDimension();
        Matrix res = copyMatrix(this);
        for (int i = 1; i < degree; i++) {
            res.multiplyByMatrix(this);
        }
        return res;
    }

    //Checks
    public boolean doesInverseExist() throws IllegalMatrixDimension {
        if(!isSquare()) return false;
        return Math.abs(this.getDeterminant()) > EPS;
    }
    public boolean isSquare(){
        return columns == rows;
    }
    public static boolean canBeMultiplied(Matrix first, Matrix second){
        return first.columns == second.rows;
    }
    public static boolean canBeAdded(Matrix first, Matrix second){
        return (first.columns == second.columns)&&
                (first.rows == second.rows);
    }

    //Auxiliary
    public static Matrix copyMatrix(Matrix matrix){
        double[][] resArray = new double[matrix.rows][matrix.columns];
        for (int i = 0; i < matrix.rows; i++) {
            if (matrix.columns >= 0) System.arraycopy(matrix.array[i], 0, resArray[i], 0, matrix.columns);
        }
        return new Matrix(resArray);
    }
    public static double[] toVector(Matrix matrix) throws IllegalMatrixDimension {
        if(matrix.rows > 1 && matrix.columns > 1) throw new IllegalMatrixDimension();
        int size = Math.max(matrix.rows, matrix.columns);

        double[] res = new double[size];
        for (int i = 0; i < size; i++) {
            if(matrix.rows == size) res[i] = matrix.array[i][0];
            else res[i] = matrix.array[0][i];
        }

        return res;
    }

    private Matrix(int rows, int columns){
        this.rows = rows;
        this.columns = columns;
    }
    private void complete(double[][] array){
        this.array = array;
    }
    private void complete(){
        this.rows = this.array.length;
        this.columns = this.array[0].length;
    }
    private void changeTo(Matrix matrix){
        this.array = matrix.array;
        this.rows = matrix.rows;
        this.columns = matrix.columns;
    }

}

class IllegalMatrixDimension extends Exception{
    /*...*/
}
